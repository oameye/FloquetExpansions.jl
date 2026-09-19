struct OpenLegPlanState{H,L}
  harmonic::H
  legs::Tuple{Vararg{L}}
end

struct OpenLegPlanVertex{H,L}
  harmonic::H
  legs::Tuple{Vararg{L}}
end

struct OpenLegResidualNode{H,L}
  order::Int
  state::OpenLegPlanState{H,L}
  generator_edges::Vector{Tuple{Int,Int,Int,OpenLegPlanState{H,L}}}
  fold_edges::Vector{Tuple{Int,OpenLegPlanState{H,L},Int,OpenLegPlanState{H,L}}}
end

struct OpenLegEvaluationPlanCounts
  residual_nodes::Int
  wave_nodes::Int
  effective_nodes::Int
  generator_products::Int
  fold_products::Int
end

struct OpenLegEvaluationPlan{H,L}
  vertices::Vector{Vector{OpenLegPlanVertex{H,L}}}
  order::Int
  zero_harmonic::H
  wave_support::Vector{Vector{OpenLegPlanState{H,L}}}
  effective_support::Vector{Vector{OpenLegPlanState{H,L}}}
  residuals::Vector{Vector{OpenLegResidualNode{H,L}}}
  counts::OpenLegEvaluationPlanCounts
end

function openleg_plan_state(harmonic::H, legs::Tuple{Vararg{L}}) where {H,L}
  return OpenLegPlanState{H,L}(harmonic, legs)
end

function openleg_combine_state(
  left_harmonic::H,
  left_legs::Tuple{Vararg{L}},
  right::OpenLegPlanState{H,L},
) where {H,L}
  return OpenLegPlanState{H,L}(left_harmonic + right.harmonic, (left_legs..., right.legs...))
end

function openleg_combine_state(
  left::OpenLegPlanState{H,L}, right::OpenLegPlanState{H,L}
) where {H,L}
  return OpenLegPlanState{H,L}(
    left.harmonic + right.harmonic, (left.legs..., right.legs...)
  )
end

function openleg_push_unique!(values::Vector{T}, value::T) where {T}
  value in values || push!(values, value)
  return values
end

function compile_openleg_evaluation_plan(
  vertices::Vector{Vector{OpenLegPlanVertex{H,L}}},
  order::Int;
  zero_harmonic::H,
) where {H,L}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  isempty(vertices) && throw(ArgumentError("at least one vertex order is required"))

  wave_support = [OpenLegPlanState{H,L}[] for _ in 1:order]
  effective_support = [OpenLegPlanState{H,L}[] for _ in 1:order]
  residuals = [OpenLegResidualNode{H,L}[] for _ in 1:order]
  identity_state = OpenLegPlanState{H,L}(zero_harmonic, ())

  residual_node_count = 0
  generator_product_count = 0
  fold_product_count = 0

  for n in 1:order
    generator_edges = Dict{
      OpenLegPlanState{H,L},Vector{Tuple{Int,Int,Int,OpenLegPlanState{H,L}}}
    }()
    fold_edges = Dict{
      OpenLegPlanState{H,L},
      Vector{Tuple{Int,OpenLegPlanState{H,L},Int,OpenLegPlanState{H,L}}},
    }()
    residual_support = OpenLegPlanState{H,L}[]

    for vertex_order in 1:min(n, length(vertices))
      previous_order = n - vertex_order
      previous_support = iszero(previous_order) ? [identity_state] : wave_support[previous_order]
      for (vertex_index, vertex) in enumerate(vertices[vertex_order]), previous in previous_support
        output = openleg_combine_state(vertex.harmonic, vertex.legs, previous)
        openleg_push_unique!(residual_support, output)
        edges = get!(generator_edges, output) do
          Tuple{Int,Int,Int,OpenLegPlanState{H,L}}[]
        end
        push!(edges, (vertex_order, vertex_index, previous_order, previous))
        generator_product_count += 1
      end
    end

    for wave_order in 1:(n - 1)
      effective_order = n - wave_order
      for wave_state in wave_support[wave_order], effective_state in effective_support[effective_order]
        output = openleg_combine_state(wave_state, effective_state)
        openleg_push_unique!(residual_support, output)
        edges = get!(fold_edges, output) do
          Tuple{Int,OpenLegPlanState{H,L},Int,OpenLegPlanState{H,L}}[]
        end
        push!(edges, (wave_order, wave_state, effective_order, effective_state))
        fold_product_count += 1
      end
    end

    for state in residual_support
      push!(
        residuals[n],
        OpenLegResidualNode(
          n,
          state,
          get(generator_edges, state, Tuple{Int,Int,Int,OpenLegPlanState{H,L}}[]),
          get(
            fold_edges,
            state,
            Tuple{Int,OpenLegPlanState{H,L},Int,OpenLegPlanState{H,L}}[],
          ),
        ),
      )
      residual_node_count += 1
    end

    effective_support[n] = OpenLegPlanState{H,L}[
      state for state in residual_support if iszero(state.harmonic)
    ]
    wave_support[n] = OpenLegPlanState{H,L}[
      state for state in residual_support if !iszero(state.harmonic)
    ]
  end

  counts = OpenLegEvaluationPlanCounts(
    residual_node_count,
    sum(length, wave_support),
    sum(length, effective_support),
    generator_product_count,
    fold_product_count,
  )
  return OpenLegEvaluationPlan(
    vertices,
    order,
    zero_harmonic,
    wave_support,
    effective_support,
    residuals,
    counts,
  )
end

function evaluate_openleg_evaluation_plan(
  plan::OpenLegEvaluationPlan{H,L},
  values::Vector{Vector{T}};
  product,
  identity_component::T,
  zero_component::T,
  inverse_weight,
) where {H,L,T}
  length(values) == length(plan.vertices) ||
    throw(ArgumentError("vertex metadata/value orders are inconsistent"))
  all(length(values[r]) == length(plan.vertices[r]) for r in eachindex(values)) ||
    throw(ArgumentError("vertex metadata/value counts are inconsistent"))

  wave = Vector{Dict{OpenLegPlanState{H,L},T}}()
  effective = Vector{Dict{OpenLegPlanState{H,L},T}}()
  identity_state = OpenLegPlanState{H,L}(plan.zero_harmonic, ())

  for n in 1:plan.order
    residual = Dict{OpenLegPlanState{H,L},T}()
    for node in plan.residuals[n]
      value = zero_component
      for (vertex_order, vertex_index, previous_order, previous_state) in node.generator_edges
        previous = if iszero(previous_order)
          previous_state == identity_state ||
            throw(ArgumentError("order-zero state must be the identity state"))
          identity_component
        else
          wave[previous_order][previous_state]
        end
        value += product(values[vertex_order][vertex_index], previous)
      end
      for (wave_order, wave_state, effective_order, effective_state) in node.fold_edges
        value -= product(wave[wave_order][wave_state], effective[effective_order][effective_state])
      end
      residual[node.state] = value
    end

    Bn = Dict{OpenLegPlanState{H,L},T}()
    Xn = Dict{OpenLegPlanState{H,L},T}()
    for (state, value) in residual
      if iszero(state.harmonic)
        Bn[state] = value
      else
        Xn[state] = inverse_weight(state.harmonic) * value
      end
    end
    push!(effective, Bn)
    push!(wave, Xn)
  end

  return (; wave, effective)
end

function collapse_openleg_states(
  components::AbstractDict{OpenLegPlanState{H,L},T}, zero_component::T
) where {H,L,T}
  collapsed = Dict{Tuple{Int,Int},T}()
  for (state, value) in components
    key = (state.harmonic, length(state.legs))
    collapsed[key] = get(collapsed, key, zero_component) + value
  end
  return openleg_periodic(collapsed, zero_component)
end

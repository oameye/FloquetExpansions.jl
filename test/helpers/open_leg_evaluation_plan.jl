struct OpenLegOutputLeg{C,S}
  channel::C
  sideband::S
end

struct OpenLegPlanState{H}
  harmonic::H
  path::Int
  grade::Int
end

struct OpenLegPlanVertex{H}
  order::Int
  harmonic::H
  path::Int
  grade::Int
  component::Int
end

struct OpenLegGeneratorEdge{H}
  vertex::Int
  wave_order::Int
  wave_state::OpenLegPlanState{H}
end

struct OpenLegFoldEdge{H}
  wave_order::Int
  wave_state::OpenLegPlanState{H}
  effective_order::Int
  effective_state::OpenLegPlanState{H}
end

struct OpenLegResidualNode{H}
  order::Int
  state::OpenLegPlanState{H}
  generator_edges::Vector{OpenLegGeneratorEdge{H}}
  fold_edges::Vector{OpenLegFoldEdge{H}}
end

struct OpenLegEvaluationPlanCounts
  residual_nodes::Int
  wave_nodes::Int
  effective_nodes::Int
  generator_products::Int
  fold_products::Int
  symbolic_paths::Int
end

struct OpenLegCompiledInput{H,T}
  vertices_by_order::Vector{Vector{OpenLegPlanVertex{H}}}
  values::Vector{T}
  paths::Vector{Tuple}
  zero_component::T
end

struct OpenLegEvaluationPlan{H,T}
  input::OpenLegCompiledInput{H,T}
  order::Int
  zero_harmonic::H
  wave_support::Vector{Vector{OpenLegPlanState{H}}}
  effective_support::Vector{Vector{OpenLegPlanState{H}}}
  residuals::Vector{Vector{OpenLegResidualNode{H}}}
  counts::OpenLegEvaluationPlanCounts
end

function openleg_intern_path!(paths::Vector{Tuple}, path_ids::Dict{Tuple,Int}, path::Tuple)
  return get!(path_ids, path) do
    push!(paths, path)
    return length(paths)
  end
end

function openleg_compile_input(
  amplitude_orders::Vector{P}; leg_factory
) where {T,P<:OpenLegPeriodic{T}}
  isempty(amplitude_orders) &&
    throw(ArgumentError("at least one amplitude order is required"))

  paths = Tuple[()]
  path_ids = Dict{Tuple,Int}(() => 1)
  values = T[]
  vertices_by_order = Vector{Vector{OpenLegPlanVertex{Int}}}()

  for (order, amplitude) in enumerate(amplitude_orders)
    vertices = OpenLegPlanVertex{Int}[]
    for ((harmonic, grade), component) in amplitude.components
      legs = leg_factory(order, harmonic, grade)
      legs isa Tuple || throw(ArgumentError("leg_factory must return a tuple"))
      length(legs) == grade ||
        throw(ArgumentError("leg_factory must return exactly grade external legs"))
      path = openleg_intern_path!(paths, path_ids, legs)
      push!(values, component)
      push!(
        vertices,
        OpenLegPlanVertex(order, harmonic, path, grade, length(values)),
      )
    end
    push!(vertices_by_order, vertices)
  end

  return OpenLegCompiledInput(
    vertices_by_order, values, paths, first(amplitude_orders).zero_component
  )
end

function openleg_compose_state!(
  left_harmonic,
  left_path::Int,
  left_grade::Int,
  right::OpenLegPlanState{H},
  paths::Vector{Tuple},
  path_ids::Dict{Tuple,Int},
) where {H}
  path = openleg_intern_path!(
    paths, path_ids, (paths[left_path]..., paths[right.path]...)
  )
  return OpenLegPlanState(left_harmonic + right.harmonic, path, left_grade + right.grade)
end

function openleg_compose_states!(
  left::OpenLegPlanState{H},
  right::OpenLegPlanState{H},
  paths::Vector{Tuple},
  path_ids::Dict{Tuple,Int},
) where {H}
  return openleg_compose_state!(
    left.harmonic, left.path, left.grade, right, paths, path_ids
  )
end

function openleg_push_unique_state!(
  states::Vector{OpenLegPlanState{H}}, state::OpenLegPlanState{H}
) where {H}
  state in states || push!(states, state)
  return states
end

openleg_is_model_state(state::OpenLegPlanState) = iszero(state.harmonic)

function compile_openleg_evaluation_plan(
  amplitude_orders::Vector{P}, order::Int; leg_factory
) where {T,P<:OpenLegPeriodic{T}}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  input = openleg_compile_input(amplitude_orders; leg_factory)
  zero_harmonic = 0
  paths = input.paths
  path_ids = Dict(path => index for (index, path) in enumerate(paths))
  identity_state = OpenLegPlanState(zero_harmonic, 1, 0)

  wave_support = [OpenLegPlanState{Int}[] for _ in 1:order]
  effective_support = [OpenLegPlanState{Int}[] for _ in 1:order]
  residuals = [OpenLegResidualNode{Int}[] for _ in 1:order]

  residual_node_count = 0
  generator_product_count = 0
  fold_product_count = 0

  for n in 1:order
    generator_edges = Dict{OpenLegPlanState{Int},Vector{OpenLegGeneratorEdge{Int}}}()
    fold_edges = Dict{OpenLegPlanState{Int},Vector{OpenLegFoldEdge{Int}}}()
    residual_support = OpenLegPlanState{Int}[]

    for vertex_order in 1:min(n, length(input.vertices_by_order))
      previous_order = n - vertex_order
      previous_states =
        iszero(previous_order) ? OpenLegPlanState{Int}[identity_state] :
        wave_support[previous_order]
      for vertex in input.vertices_by_order[vertex_order], previous in previous_states
        output = openleg_compose_state!(
          vertex.harmonic,
          vertex.path,
          vertex.grade,
          previous,
          paths,
          path_ids,
        )
        openleg_push_unique_state!(residual_support, output)
        edges = get!(generator_edges, output) do
          return OpenLegGeneratorEdge{Int}[]
        end
        push!(edges, OpenLegGeneratorEdge(vertex.component, previous_order, previous))
        generator_product_count += 1
      end
    end

    for wave_order in 1:(n - 1)
      effective_order = n - wave_order
      for wave_state in wave_support[wave_order],
        effective_state in effective_support[effective_order]

        output = openleg_compose_states!(
          wave_state, effective_state, paths, path_ids
        )
        openleg_push_unique_state!(residual_support, output)
        edges = get!(fold_edges, output) do
          return OpenLegFoldEdge{Int}[]
        end
        push!(
          edges,
          OpenLegFoldEdge(
            wave_order, wave_state, effective_order, effective_state
          ),
        )
        fold_product_count += 1
      end
    end

    for state in residual_support
      push!(
        residuals[n],
        OpenLegResidualNode(
          n,
          state,
          get(generator_edges, state, OpenLegGeneratorEdge{Int}[]),
          get(fold_edges, state, OpenLegFoldEdge{Int}[]),
        ),
      )
      residual_node_count += 1
    end

    effective_support[n] = OpenLegPlanState{Int}[
      state for state in residual_support if openleg_is_model_state(state)
    ]
    wave_support[n] = OpenLegPlanState{Int}[
      state for state in residual_support if !openleg_is_model_state(state)
    ]
  end

  compiled_input = OpenLegCompiledInput(
    input.vertices_by_order, input.values, paths, input.zero_component
  )
  counts = OpenLegEvaluationPlanCounts(
    residual_node_count,
    sum(length, wave_support),
    sum(length, effective_support),
    generator_product_count,
    fold_product_count,
    length(paths),
  )
  return OpenLegEvaluationPlan(
    compiled_input,
    order,
    zero_harmonic,
    wave_support,
    effective_support,
    residuals,
    counts,
  )
end

function evaluate_openleg_evaluation_plan(
  plan::OpenLegEvaluationPlan{H,T}; product, inverse_weight, identity_component::T
) where {H,T}
  wave = Vector{Dict{OpenLegPlanState{H},T}}()
  effective = Vector{Dict{OpenLegPlanState{H},T}}()

  for n in 1:plan.order
    residual = Dict{OpenLegPlanState{H},T}()
    for node in plan.residuals[n]
      value = plan.input.zero_component
      for edge in node.generator_edges
        left = plan.input.values[edge.vertex]
        right =
          iszero(edge.wave_order) ? identity_component :
          wave[edge.wave_order][edge.wave_state]
        value += product(left, right)
      end
      for edge in node.fold_edges
        value -= product(
          wave[edge.wave_order][edge.wave_state],
          effective[edge.effective_order][edge.effective_state],
        )
      end
      residual[node.state] = value
    end

    Bn = Dict{OpenLegPlanState{H},T}()
    Xn = Dict{OpenLegPlanState{H},T}()
    for (state, value) in residual
      if openleg_is_model_state(state)
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

function collapse_openleg_plan_component(
  components::Dict{OpenLegPlanState{H},T},
  plan::OpenLegEvaluationPlan{H,T},
) where {H,T}
  collapsed = Dict{Tuple{Int,Int},T}()
  for (state, value) in components
    key = (state.harmonic, state.grade)
    collapsed[key] = get(collapsed, key, plan.input.zero_component) + value
  end
  return openleg_periodic(collapsed, plan.input.zero_component)
end

function openleg_plan_path(plan::OpenLegEvaluationPlan, state::OpenLegPlanState)
  return plan.input.paths[state.path]
end

function openleg_first_return_mismatch(states)
  isempty(states) && return false
  total = zero(first(states).harmonic)
  for (index, state) in enumerate(states)
    total += state.harmonic
    index < length(states) && iszero(total) && return false
  end
  return iszero(total)
end

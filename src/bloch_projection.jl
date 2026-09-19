struct BlochResidualPlanNode{H}
  order::Int
  harmonic::H
  generator_edges::Vector{Tuple{H,H}}
  fold_edges::Vector{Tuple{Int,Int}}
end

struct BlochProjectionPlanCounts
  residual_nodes::Int
  wave_nodes::Int
  effective_nodes::Int
  generator_products::Int
  fold_products::Int
end

struct BlochProjectionPlan{H}
  input_support::Vector{H}
  order::Int
  zero_harmonic::H
  wave_support::Vector{Vector{H}}
  effective_active::Vector{Bool}
  residuals::Vector{Vector{BlochResidualPlanNode{H}}}
  counts::BlochProjectionPlanCounts
end

struct BlochProjectionResult{H,T}
  wave::Vector{Dict{H,T}}
  effective::Vector{T}
end

struct BlochOrderOperations{P,M,Q,S}
  product::P
  project_model::M
  solve_complement::Q
  simplifier::S
end

function BlochOrderOperations(product, project_model, solve_complement; simplifier=identity)
  return BlochOrderOperations(product, project_model, solve_complement, simplifier)
end

struct BlochOrderRecurrenceResult{T}
  wave::Vector{T}
  effective::Vector{T}
  products::Int
end

function evaluate_bloch_order_recurrence(
  generators_by_order::AbstractVector{T},
  order::Int,
  identity_component::T,
  zero_component::T,
  operations::BlochOrderOperations,
) where {T}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  isempty(generators_by_order) &&
    throw(ArgumentError("at least one generator coefficient is required"))

  wave = Vector{T}(undef, max(order - 1, 0))
  effective = Vector{T}(undef, order)
  products = 0

  for n in 1:order
    residual = zero_component

    for generator_order in 1:min(n, length(generators_by_order))
      previous_order = n - generator_order
      previous = iszero(previous_order) ? identity_component : wave[previous_order]
      residual += operations.product(generators_by_order[generator_order], previous)
      products += 1
    end

    for wave_order in 1:(n - 1)
      effective_order = n - wave_order
      residual -= operations.product(wave[wave_order], effective[effective_order])
      products += 1
    end

    residual = operations.simplifier(residual)::T
    effective[n] = operations.simplifier(operations.project_model(residual))::T
    if n < order
      wave[n] = operations.simplifier(operations.solve_complement(residual))::T
    end
  end

  return BlochOrderRecurrenceResult(wave, effective, products)
end

function bloch_push_unique!(values::Vector{H}, value::H) where {H}
  value in values || push!(values, value)
  return values
end

function compile_bloch_projection_plan(support_input, order::Int)
  support = unique(collect(support_input))
  isempty(support) && throw(ArgumentError("harmonic support must not be empty"))
  return compile_bloch_projection_plan(support, order, zero(first(support)))
end

function compile_bloch_projection_plan(
  support_input, order::Int, zero_harmonic::H
) where {H}
  order >= 1 || throw(ArgumentError("order must be >= 1"))

  support = unique(collect(support_input))
  isempty(support) && throw(ArgumentError("harmonic support must not be empty"))
  all(harmonic isa H for harmonic in support) ||
    throw(ArgumentError("harmonic support and zero harmonic must have one common type"))
  typed_support = H[harmonic for harmonic in support]

  wave_support = [H[] for _ in 1:max(order - 1, 0)]
  effective_active = falses(order)
  residuals = [BlochResidualPlanNode{H}[] for _ in 1:max(order - 1, 0)]

  effective_active[1] = zero_harmonic in typed_support
  if order > 1
    wave_support[1] = H[harmonic for harmonic in typed_support if !iszero(harmonic)]
  end

  residual_node_count = 0
  generator_product_count = 0
  fold_product_count = 0

  for n in 1:(order - 1)
    generator_edges = Dict{H,Vector{Tuple{H,H}}}()
    fold_edges = Dict{H,Vector{Tuple{Int,Int}}}()
    residual_support = H[]

    for generator_harmonic in typed_support, wave_harmonic in wave_support[n]
      output_harmonic = generator_harmonic + wave_harmonic
      bloch_push_unique!(residual_support, output_harmonic)
      edges = get!(generator_edges, output_harmonic) do
        return Tuple{H,H}[]
      end
      push!(edges, (generator_harmonic, wave_harmonic))
      generator_product_count += 1
    end

    for wave_order in 1:n
      effective_order = n - wave_order
      effective_active[effective_order + 1] || continue
      for harmonic in wave_support[wave_order]
        bloch_push_unique!(residual_support, harmonic)
        edges = get!(fold_edges, harmonic) do
          return Tuple{Int,Int}[]
        end
        push!(edges, (wave_order, effective_order))
        fold_product_count += 1
      end
    end

    for harmonic in residual_support
      push!(
        residuals[n],
        BlochResidualPlanNode(
          n,
          harmonic,
          get(generator_edges, harmonic, Tuple{H,H}[]),
          get(fold_edges, harmonic, Tuple{Int,Int}[]),
        ),
      )
      residual_node_count += 1
    end

    effective_active[n + 1] = zero_harmonic in residual_support
    if n < order - 1
      wave_support[n + 1] = H[
        harmonic for harmonic in residual_support if !iszero(harmonic)
      ]
    end
  end

  counts = BlochProjectionPlanCounts(
    residual_node_count,
    sum(length, wave_support),
    count(identity, effective_active),
    generator_product_count,
    fold_product_count,
  )
  return BlochProjectionPlan(
    typed_support,
    order,
    zero_harmonic,
    wave_support,
    collect(effective_active),
    residuals,
    counts,
  )
end

function evaluate_bloch_projection_plan(
  plan::BlochProjectionPlan{H},
  components::AbstractDict{H,T};
  product,
  inverse_weight,
  zero_component::T,
  simplifier=identity,
) where {H,T}
  wave = Vector{Dict{H,T}}()
  effective = T[]

  B0 = simplifier(get(components, plan.zero_harmonic, zero_component))::T
  push!(effective, B0)

  if plan.order > 1
    X1 = Dict{H,T}()
    for harmonic in plan.wave_support[1]
      X1[harmonic] = simplifier(inverse_weight(harmonic) * components[harmonic])::T
    end
    push!(wave, X1)
  end

  for n in 1:(plan.order - 1)
    residual = Dict{H,T}()
    for node in plan.residuals[n]
      value = zero_component
      for (generator_harmonic, wave_harmonic) in node.generator_edges
        value += product(components[generator_harmonic], wave[n][wave_harmonic])
      end
      for (wave_order, effective_order) in node.fold_edges
        value -= product(wave[wave_order][node.harmonic], effective[effective_order + 1])
      end
      residual[node.harmonic] = simplifier(value)::T
    end

    Bn = simplifier(get(residual, plan.zero_harmonic, zero_component))::T
    push!(effective, Bn)

    if n < plan.order - 1
      Xnext = Dict{H,T}()
      for harmonic in plan.wave_support[n + 1]
        Xnext[harmonic] = simplifier(inverse_weight(harmonic) * residual[harmonic])::T
      end
      push!(wave, Xnext)
    end
  end

  return BlochProjectionResult(wave, effective)
end

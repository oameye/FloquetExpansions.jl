struct BlochResidualNode{H}
  order::Int
  harmonic::H
  generator_edges::Vector{Tuple{H,H}}
  fold_edges::Vector{Tuple{Int,Int}}
end

struct BlochEvaluationPlanCounts
  residual_nodes::Int
  wave_nodes::Int
  effective_nodes::Int
  generator_products::Int
  fold_products::Int
end

struct BlochEvaluationPlan{H}
  input_support::Vector{H}
  order::Int
  zero_harmonic::H
  wave_support::Vector{Vector{H}}
  effective_active::Vector{Bool}
  residuals::Vector{Vector{BlochResidualNode{H}}}
  counts::BlochEvaluationPlanCounts
end

function push_unique!(values::Vector{H}, value::H) where {H}
  value in values || push!(values, value)
  return values
end

function compile_bloch_evaluation_plan(
  support_input, order::Int; zero_harmonic=zero(first(support_input))
)
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  isempty(support_input) && throw(ArgumentError("harmonic support must not be empty"))

  support = unique(collect(support_input))
  H = typeof(zero_harmonic)
  all(harmonic isa H for harmonic in support) ||
    throw(ArgumentError("harmonic support and zero harmonic must have one common type"))
  support = H[harmonic for harmonic in support]

  wave_support = [H[] for _ in 1:max(order - 1, 0)]
  effective_active = falses(order)
  residuals = [BlochResidualNode{H}[] for _ in 1:max(order - 1, 0)]

  effective_active[1] = zero_harmonic in support
  if order > 1
    wave_support[1] = H[harmonic for harmonic in support if !iszero(harmonic)]
  end

  residual_node_count = 0
  generator_product_count = 0
  fold_product_count = 0

  for n in 1:(order - 1)
    generator_edges = Dict{H,Vector{Tuple{H,H}}}()
    fold_edges = Dict{H,Vector{Tuple{Int,Int}}}()
    residual_support = H[]

    for generator_harmonic in support, wave_harmonic in wave_support[n]
      output_harmonic = generator_harmonic + wave_harmonic
      push_unique!(residual_support, output_harmonic)
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
        push_unique!(residual_support, harmonic)
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
        BlochResidualNode(
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

  counts = BlochEvaluationPlanCounts(
    residual_node_count,
    sum(length, wave_support),
    count(identity, effective_active),
    generator_product_count,
    fold_product_count,
  )
  return BlochEvaluationPlan(
    support,
    order,
    zero_harmonic,
    wave_support,
    collect(effective_active),
    residuals,
    counts,
  )
end

function evaluate_bloch_evaluation_plan(
  plan::BlochEvaluationPlan{H},
  components::AbstractDict{H,T};
  product,
  inverse_weight,
  zero_component::T,
  simplifier=identity,
) where {H,T}
  wave = Vector{Dict{H,T}}()
  effective = T[]

  B0 = simplifier(get(components, plan.zero_harmonic, zero_component))
  push!(effective, B0)

  if plan.order > 1
    X1 = Dict{H,T}()
    for harmonic in plan.wave_support[1]
      X1[harmonic] = simplifier(inverse_weight(harmonic) * components[harmonic])
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
      residual[node.harmonic] = simplifier(value)
    end

    Bn = simplifier(get(residual, plan.zero_harmonic, zero_component))
    push!(effective, Bn)

    if n < plan.order - 1
      Xnext = Dict{H,T}()
      for harmonic in plan.wave_support[n + 1]
        Xnext[harmonic] = simplifier(inverse_weight(harmonic) * residual[harmonic])
      end
      push!(wave, Xnext)
    end
  end

  return (; wave, effective)
end

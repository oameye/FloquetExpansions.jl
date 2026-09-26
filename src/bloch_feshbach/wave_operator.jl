struct BlochConventions{F,W<:Number,K<:Number}
  product::F
  weight_phase::W
  kick_phase::K
end

struct WaveOperatorNode{H}
  harmonic::H
  generator_edges::Vector{Tuple{H,H}}
  fold_edges::Vector{Tuple{Int,Int}}
end

struct WaveOperatorPlan{H}
  support::Vector{H}
  order::Int
  zero_harmonic::H
  wave_support::Vector{Vector{H}}
  residuals::Vector{Vector{WaveOperatorNode{H}}}
end

struct WaveOperator{H,T}
  coefficients::Vector{Dict{H,T}}
  bloch_effective_generator::Vector{T}
end

bloch_harmonic_inverse(harmonic::Int) = 1 // harmonic

simplify_component(component) = SQA.simplify(component)

function push_unique!(values::Vector{H}, value::H) where {H}
  value in values || push!(values, value)
  return values
end

function nonzero_harmonics(harmonics::Vector{H}, zero_harmonic::H) where {H}
  return H[harmonic for harmonic in harmonics if harmonic != zero_harmonic]
end

function compile_residual_nodes(
  support::Vector{H},
  wave_support::Vector{Vector{H}},
  bloch_effective_active::BitVector,
  n::Int,
) where {H}
  generator_edges = Dict{H,Vector{Tuple{H,H}}}()
  fold_edges = Dict{H,Vector{Tuple{Int,Int}}}()
  residual_support = H[]

  for generator_harmonic in support, wave_harmonic in wave_support[n]
    harmonic = generator_harmonic + wave_harmonic
    push_unique!(residual_support, harmonic)
    edges = get!(Vector{Tuple{H,H}}, generator_edges, harmonic)
    push!(edges, (generator_harmonic, wave_harmonic))
  end

  for wave_order in 1:n
    effective_order = n - wave_order
    bloch_effective_active[effective_order + 1] || continue
    for harmonic in wave_support[wave_order]
      push_unique!(residual_support, harmonic)
      edges = get!(Vector{Tuple{Int,Int}}, fold_edges, harmonic)
      push!(edges, (wave_order, effective_order))
    end
  end

  return WaveOperatorNode{H}[
    WaveOperatorNode(
      harmonic,
      get(generator_edges, harmonic, Tuple{H,H}[]),
      get(fold_edges, harmonic, Tuple{Int,Int}[]),
    ) for harmonic in residual_support
  ]
end

function compile_wave_operator_plan(support_input, order::Int, zero_harmonic::H) where {H}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  support = unique(H[harmonic for harmonic in support_input])
  isempty(support) && throw(ArgumentError("harmonic support must not be empty"))

  wave_support = [H[] for _ in 1:(order - 1)]
  residuals = [WaveOperatorNode{H}[] for _ in 1:(order - 1)]
  bloch_effective_active = falses(order)
  bloch_effective_active[1] = zero_harmonic in support
  order > 1 && (wave_support[1] = nonzero_harmonics(support, zero_harmonic))

  for n in 1:(order - 1)
    residuals[n] = compile_residual_nodes(support, wave_support, bloch_effective_active, n)
    residual_support = H[node.harmonic for node in residuals[n]]
    bloch_effective_active[n + 1] = zero_harmonic in residual_support
    if n < order - 1
      wave_support[n + 1] = nonzero_harmonics(residual_support, zero_harmonic)
    end
  end

  return WaveOperatorPlan(support, order, zero_harmonic, wave_support, residuals)
end

function bloch_antiderivative(
  values::AbstractDict{H,T}, harmonics::Vector{H}, weight_phase
) where {H,T}
  result = Dict{H,T}()
  for harmonic in harmonics
    weight = weight_phase * bloch_harmonic_inverse(harmonic)
    result[harmonic] = simplify_component(weight * values[harmonic])::T
  end
  return result
end

function evaluate_residual(
  node::WaveOperatorNode{H},
  components::AbstractDict{H,T},
  coefficients::Vector{Dict{H,T}},
  bloch_effective_generator::Vector{T},
  n::Int,
  zero_component::T,
  product,
) where {H,T}
  value = zero_component
  for (generator_harmonic, wave_harmonic) in node.generator_edges
    value += product(components[generator_harmonic], coefficients[n][wave_harmonic])
  end
  for (wave_order, effective_order) in node.fold_edges
    wave = coefficients[wave_order][node.harmonic]
    value -= product(wave, bloch_effective_generator[effective_order + 1])
  end
  return simplify_component(value)::T
end

function evaluate_wave_operator(
  plan::WaveOperatorPlan{H},
  components::AbstractDict{H,T},
  zero_component::T,
  conventions::BlochConventions,
) where {H,T}
  (; product, weight_phase) = conventions
  coefficients = Dict{H,T}[]
  static_part = get(components, plan.zero_harmonic, zero_component)
  bloch_effective_generator = T[simplify_component(static_part)::T]
  if plan.order > 1
    push!(
      coefficients, bloch_antiderivative(components, plan.wave_support[1], weight_phase)
    )
  end

  for n in 1:(plan.order - 1)
    residual = Dict{H,T}()
    for node in plan.residuals[n]
      residual[node.harmonic] = evaluate_residual(
        node,
        components,
        coefficients,
        bloch_effective_generator,
        n,
        zero_component,
        product,
      )
    end
    push!(bloch_effective_generator, get(residual, plan.zero_harmonic, zero_component))
    if n < plan.order - 1
      push!(
        coefficients, bloch_antiderivative(residual, plan.wave_support[n + 1], weight_phase)
      )
    end
  end

  return WaveOperator(coefficients, bloch_effective_generator)
end

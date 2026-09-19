if !isdefined(@__MODULE__, :MixedQRVertex)
  include(joinpath(@__DIR__, "mixed_qr_kernel_reference.jl"))
end

struct GramQRVertex{T}
  metadata::MixedQRVertex
  value::T
end

struct GramStinespringWord{T}
  metadata::Vector{MixedQRVertex}
  output_number::Int
  drift_order::Int
  system_value::T
end

struct FormalPeriodMatrix{T}
  coefficients::Vector{T}
end

function formal_period_zero(zero_component::T) where {T}
  return FormalPeriodMatrix{T}(T[zero_component])
end

function formal_period_accumulate!(
  polynomial::FormalPeriodMatrix{T},
  scalar_polynomial::Vector{S},
  matrix::T,
) where {S,T}
  required = length(scalar_polynomial)
  while length(polynomial.coefficients) < required
    push!(polynomial.coefficients, zero(matrix))
  end
  for index in eachindex(scalar_polynomial)
    polynomial.coefficients[index] += scalar_polynomial[index] * matrix
  end
  return polynomial
end

function formal_period_coefficient(
  polynomial::FormalPeriodMatrix{T}, degree::Int, zero_component::T
) where {T}
  degree >= 0 || throw(ArgumentError("degree must be nonnegative"))
  degree + 1 <= length(polynomial.coefficients) || return zero_component
  return polynomial.coefficients[degree + 1]
end

function gram_qr_jump(channel::Symbol, harmonic::Int, value::T) where {T}
  return GramQRVertex(mixed_qr_jump(channel, harmonic), value)
end

function gram_qr_drift(harmonic::Int, value::T) where {T}
  return GramQRVertex(mixed_qr_drift(harmonic), value)
end

function gram_stinespring_triangle_count(
  jump_label_count::Int, drift_label_count::Int, generator_order::Int
)
  jump_label_count >= 0 || throw(ArgumentError("jump label count must be nonnegative"))
  drift_label_count >= 0 || throw(ArgumentError("drift label count must be nonnegative"))
  generator_order >= 0 || throw(ArgumentError("generator order must be nonnegative"))
  alphabet = jump_label_count + drift_label_count
  return sum(alphabet^length for length in 0:(generator_order + 1))
end

function gram_stinespring_qd_count(
  jump_label_count::Int,
  drift_label_count::Int,
  output_number::Int,
  drift_order::Int,
)
  output_number >= 0 || throw(ArgumentError("output number must be nonnegative"))
  drift_order >= 0 || throw(ArgumentError("drift order must be nonnegative"))
  vertex_count = output_number + drift_order
  return binomial(vertex_count, output_number) *
         jump_label_count^output_number *
         drift_label_count^drift_order
end

function gram_stinespring_triangle_words(
  jumps::Vector{GramQRVertex{T}},
  drifts::Vector{GramQRVertex{T}},
  generator_order::Int,
  identity_component::T,
) where {T}
  generator_order >= 0 || throw(ArgumentError("generator order must be nonnegative"))
  all(mixed_qr_is_jump(vertex.metadata) for vertex in jumps) ||
    throw(ArgumentError("jump alphabet contains a drift vertex"))
  all(!mixed_qr_is_jump(vertex.metadata) for vertex in drifts) ||
    throw(ArgumentError("drift alphabet contains a jump vertex"))

  alphabet = GramQRVertex{T}[jumps; drifts]
  vacuum = GramStinespringWord(MixedQRVertex[], 0, 0, identity_component)
  words = GramStinespringWord{T}[vacuum]
  frontier = GramStinespringWord{T}[vacuum]

  for _ in 1:(generator_order + 1)
    next_frontier = GramStinespringWord{T}[]
    for prefix in frontier, vertex in alphabet
      metadata = [prefix.metadata; vertex.metadata]
      output_number = prefix.output_number + Int(mixed_qr_is_jump(vertex.metadata))
      drift_order = prefix.drift_order + Int(!mixed_qr_is_jump(vertex.metadata))
      system_value = vertex.value * prefix.system_value
      push!(
        next_frontier,
        GramStinespringWord(metadata, output_number, drift_order, system_value),
      )
    end
    append!(words, next_frontier)
    frontier = next_frontier
  end
  return words
end

function gram_stinespring_metric_series(
  words::Vector{GramStinespringWord{T}}, generator_order::Int, zero_component::T
) where {T}
  generator_order >= 0 || throw(ArgumentError("generator order must be nonnegative"))
  max_channel_order = 2 * (generator_order + 1)
  metric = [formal_period_zero(zero_component) for _ in 0:max_channel_order]

  for left in words, right in words
    left.output_number == right.output_number || continue
    channel_order = left.output_number + left.drift_order + right.drift_order
    gram = mixed_qr_gram_polynomial(left.metadata, right.metadata)
    all(iszero, gram) && continue
    system_term = adjoint(left.system_value) * right.system_value
    formal_period_accumulate!(metric[channel_order + 1], gram, system_term)
  end
  return metric
end

function gram_stinespring_channel_series(
  words::Vector{GramStinespringWord{T}},
  generator_order::Int,
  zero_superoperator::S;
  paste,
) where {T,S}
  max_channel_order = 2 * (generator_order + 1)
  channel = [formal_period_zero(zero_superoperator) for _ in 0:max_channel_order]

  for left in words, right in words
    left.output_number == right.output_number || continue
    channel_order = left.output_number + left.drift_order + right.drift_order
    gram = mixed_qr_gram_polynomial(left.metadata, right.metadata)
    all(iszero, gram) && continue
    system_term = paste(left.system_value, right.system_value)
    formal_period_accumulate!(channel[channel_order + 1], conj.(gram), system_term)
  end
  return channel
end

function formal_period_iszero(polynomial::FormalPeriodMatrix)
  return all(iszero, polynomial.coefficients)
end

function formal_period_isapprox_zero(polynomial::FormalPeriodMatrix; atol)
  return all(coefficient -> all(abs.(coefficient) .<= atol), polynomial.coefficients)
end

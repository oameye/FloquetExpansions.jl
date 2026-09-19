struct CKEndpointConstraint
  start::Int
  stop::Int

  function CKEndpointConstraint(start::Int, stop::Int)
    0 <= start < stop ||
      throw(ArgumentError("endpoint constraint must span a nonempty interval"))
    return new(start, stop)
  end
end

struct CKEndpointKey{H}
  vertices::Vector{CKBranchVertex{H}}
  model_intervals::Vector{CKEndpointConstraint}
  resolvent_intervals::Vector{CKEndpointConstraint}
end

function Base.:(==)(left::CKEndpointKey, right::CKEndpointKey)
  return left.vertices == right.vertices &&
         left.model_intervals == right.model_intervals &&
         left.resolvent_intervals == right.resolvent_intervals
end

Base.isequal(left::CKEndpointKey, right::CKEndpointKey) = left == right

function Base.hash(key::CKEndpointKey, seed::UInt)
  result = hash(length(key.vertices), seed)
  for vertex in key.vertices
    result = hash(vertex, result)
  end
  result = hash(length(key.model_intervals), result)
  for interval in key.model_intervals
    result = hash(interval, result)
  end
  result = hash(length(key.resolvent_intervals), result)
  for interval in key.resolvent_intervals
    result = hash(interval, result)
  end
  return result
end

struct CKEndpointKernel{H,T}
  terms::Dict{CKEndpointKey{H},T}
  zero_component::T
end

function ck_endpoint_accumulate!(
  terms::Dict{CKEndpointKey{H},T}, key::CKEndpointKey{H}, value::T, zero_component::T
) where {H,T}
  updated = get(terms, key, zero_component) + value
  if ck_kernel_iszero(updated)
    haskey(terms, key) && delete!(terms, key)
  else
    terms[key] = updated
  end
  return terms
end

function ck_endpoint_kernel(state::CKPhysicalKernel{H,T}) where {H,T}
  result = Dict{CKEndpointKey{H},T}()
  for (key, value) in state.terms
    model_intervals = CKEndpointConstraint[
      CKEndpointConstraint(0, cut) for cut in key.model_cuts
    ]
    resolvent_intervals = CKEndpointConstraint[
      CKEndpointConstraint(0, cut) for cut in key.resolvent_cuts
    ]
    endpoint_key = CKEndpointKey(copy(key.vertices), model_intervals, resolvent_intervals)
    ck_endpoint_accumulate!(result, endpoint_key, value, state.zero_component)
  end
  return CKEndpointKernel(result, state.zero_component)
end

function ck_endpoint_zero(state::CKEndpointKernel{H,T}) where {H,T}
  return CKEndpointKernel(Dict{CKEndpointKey{H},T}(), state.zero_component)
end

function ck_endpoint_shifted_intervals(intervals::Vector{CKEndpointConstraint}, shift::Int)
  return CKEndpointConstraint[
    CKEndpointConstraint(interval.start + shift, interval.stop + shift) for
    interval in intervals
  ]
end

function ck_endpoint_product(
  left::CKEndpointKernel{H,T}, right::CKEndpointKernel{H,T}
) where {H,T}
  result = Dict{CKEndpointKey{H},T}()
  for (left_key, left_value) in left.terms, (right_key, right_value) in right.terms
    right_length = length(right_key.vertices)
    vertices = CKBranchVertex{H}[right_key.vertices; left_key.vertices]
    model_intervals = CKEndpointConstraint[
      right_key.model_intervals;
      ck_endpoint_shifted_intervals(left_key.model_intervals, right_length)
    ]
    resolvent_intervals = CKEndpointConstraint[
      right_key.resolvent_intervals;
      ck_endpoint_shifted_intervals(left_key.resolvent_intervals, right_length)
    ]
    key = CKEndpointKey(vertices, model_intervals, resolvent_intervals)
    ck_endpoint_accumulate!(result, key, left_value * right_value, left.zero_component)
  end
  return CKEndpointKernel(result, left.zero_component)
end

function Base.:(==)(left::CKEndpointKernel, right::CKEndpointKernel)
  return left.zero_component == right.zero_component && left.terms == right.terms
end

function Base.:+(left::CKEndpointKernel{H,T}, right::CKEndpointKernel{H,T}) where {H,T}
  result = copy(left.terms)
  for (key, value) in right.terms
    ck_endpoint_accumulate!(result, key, value, left.zero_component)
  end
  return CKEndpointKernel(result, left.zero_component)
end

function Base.:-(left::CKEndpointKernel{H,T}, right::CKEndpointKernel{H,T}) where {H,T}
  result = copy(left.terms)
  for (key, value) in right.terms
    ck_endpoint_accumulate!(result, key, -value, left.zero_component)
  end
  return CKEndpointKernel(result, left.zero_component)
end

function ck_endpoint_divide(state::CKEndpointKernel{H,T}, divisor::Int) where {H,T}
  iszero(divisor) && throw(DivideError())
  result = Dict{CKEndpointKey{H},T}()
  for (key, value) in state.terms
    ck_endpoint_accumulate!(result, key, value / divisor, state.zero_component)
  end
  return CKEndpointKernel(result, state.zero_component)
end

struct CKPeriodPolynomial{H,T}
  coefficients::Vector{CKEndpointKernel{H,T}}
end

function ck_period_polynomial(coefficients::Vector{CKEndpointKernel{H,T}}) where {H,T}
  isempty(coefficients) &&
    throw(ArgumentError("period polynomial must contain a coefficient"))
  last_nonzero = findlast(coefficient -> !isempty(coefficient.terms), coefficients)
  length_to_keep = isnothing(last_nonzero) ? 1 : last_nonzero
  return CKPeriodPolynomial{H,T}(copy(coefficients[1:length_to_keep]))
end

function ck_period_zero(state::CKEndpointKernel{H,T}) where {H,T}
  return ck_period_polynomial(CKEndpointKernel{H,T}[ck_endpoint_zero(state)])
end

function ck_period_constant(state::CKEndpointKernel{H,T}) where {H,T}
  return ck_period_polynomial(CKEndpointKernel{H,T}[state])
end

function ck_period_monomial(state::CKEndpointKernel{H,T}, power::Int) where {H,T}
  power >= 0 || throw(ArgumentError("period power must be nonnegative"))
  coefficients = CKEndpointKernel{H,T}[ck_endpoint_zero(state) for _ in 0:power]
  coefficients[power + 1] = state
  return ck_period_polynomial(coefficients)
end

function Base.:(==)(left::CKPeriodPolynomial, right::CKPeriodPolynomial)
  return left.coefficients == right.coefficients
end

function Base.:+(left::CKPeriodPolynomial{H,T}, right::CKPeriodPolynomial{H,T}) where {H,T}
  zero_state = ck_endpoint_zero(first(left.coefficients))
  length_result = max(length(left.coefficients), length(right.coefficients))
  coefficients = Vector{CKEndpointKernel{H,T}}(undef, length_result)
  for index in 1:length_result
    left_value = index <= length(left.coefficients) ? left.coefficients[index] : zero_state
    right_value =
      index <= length(right.coefficients) ? right.coefficients[index] : zero_state
    coefficients[index] = left_value + right_value
  end
  return ck_period_polynomial(coefficients)
end

function Base.:-(left::CKPeriodPolynomial{H,T}, right::CKPeriodPolynomial{H,T}) where {H,T}
  zero_state = ck_endpoint_zero(first(left.coefficients))
  length_result = max(length(left.coefficients), length(right.coefficients))
  coefficients = Vector{CKEndpointKernel{H,T}}(undef, length_result)
  for index in 1:length_result
    left_value = index <= length(left.coefficients) ? left.coefficients[index] : zero_state
    right_value =
      index <= length(right.coefficients) ? right.coefficients[index] : zero_state
    coefficients[index] = left_value - right_value
  end
  return ck_period_polynomial(coefficients)
end

function Base.:-(value::CKPeriodPolynomial{H,T}) where {H,T}
  return ck_period_zero(first(value.coefficients)) - value
end

function Base.:*(left::CKPeriodPolynomial{H,T}, right::CKPeriodPolynomial{H,T}) where {H,T}
  zero_state = ck_endpoint_zero(first(left.coefficients))
  maximum_power = length(left.coefficients) + length(right.coefficients) - 2
  coefficients = CKEndpointKernel{H,T}[
    ck_endpoint_zero(zero_state) for _ in 0:maximum_power
  ]
  for left_index in eachindex(left.coefficients),
    right_index in eachindex(right.coefficients)

    output_index = left_index + right_index - 1
    coefficients[output_index] += ck_endpoint_product(
      left.coefficients[left_index], right.coefficients[right_index]
    )
  end
  return ck_period_polynomial(coefficients)
end

function Base.:/(value::CKPeriodPolynomial{H,T}, divisor::Int) where {H,T}
  return ck_period_polynomial([
    ck_endpoint_divide(coefficient, divisor) for coefficient in value.coefficients
  ])
end

function ck_truncated_series_product(
  left::AbstractVector{T}, right::AbstractVector{T}, order::Int, zero_value::T
) where {T}
  order >= 0 || throw(ArgumentError("series order must be nonnegative"))
  result = Vector{T}(undef, order + 1)
  for n in 0:order
    coefficient = zero_value
    for left_order in 0:n
      right_order = n - left_order
      left_order + 1 <= length(left) || continue
      right_order + 1 <= length(right) || continue
      coefficient = coefficient + left[left_order + 1] * right[right_order + 1]
    end
    result[n + 1] = coefficient
  end
  return result
end

function ck_unit_series_inverse(
  series::AbstractVector{T}, order::Int, one_value::T, zero_value::T
) where {T}
  order >= 0 || throw(ArgumentError("series order must be nonnegative"))
  isempty(series) && throw(ArgumentError("series must contain its order-zero coefficient"))
  series[1] == one_value ||
    throw(ArgumentError("series must have unit order-zero coefficient"))

  result = Vector{T}(undef, order + 1)
  result[1] = one_value
  for n in 1:order
    coefficient = zero_value
    for left_order in 1:n
      left_order + 1 <= length(series) || continue
      coefficient = coefficient + series[left_order + 1] * result[n - left_order + 1]
    end
    result[n + 1] = -coefficient
  end
  return result
end

function ck_zero_constant_series_exponential(
  series::AbstractVector{T}, order::Int, one_value::T, zero_value::T
) where {T}
  order >= 0 || throw(ArgumentError("series order must be nonnegative"))
  isempty(series) && throw(ArgumentError("series must contain its order-zero coefficient"))
  series[1] == zero_value ||
    throw(ArgumentError("exponential input must have zero order-zero coefficient"))

  result = T[zero_value for _ in 0:order]
  result[1] = one_value
  term = T[zero_value for _ in 0:order]
  term[1] = one_value

  for power in 1:order
    term = ck_truncated_series_product(term, series, order, zero_value)
    for n in 0:order
      term[n + 1] = term[n + 1] / power
      result[n + 1] = result[n + 1] + term[n + 1]
    end
  end
  return result
end

struct CKPeriodAmplitudeResult{P}
  amplitude::Vector{P}
  endpoint_wave::Vector{P}
  endpoint_wave_inverse::Vector{P}
  slow_propagator::Vector{P}
end

function evaluate_ck_period_amplitude(
  effective::AbstractVector{K},
  wave::AbstractVector{K},
  order::Int,
  identity_state::K,
  zero_state::K,
) where {K<:CKPhysicalKernel}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  length(effective) >= order ||
    throw(ArgumentError("effective series does not contain the requested order"))
  length(wave) >= order - 1 ||
    throw(ArgumentError("wave series does not contain the required endpoint orders"))

  zero_endpoint = ck_endpoint_kernel(zero_state)
  identity_endpoint = ck_endpoint_kernel(identity_state)
  zero_period = ck_period_zero(zero_endpoint)
  identity_period = ck_period_constant(identity_endpoint)

  endpoint_wave = typeof(identity_period)[zero_period for _ in 0:order]
  endpoint_wave[1] = identity_period
  for n in 1:(order - 1)
    endpoint_wave[n + 1] = ck_period_constant(ck_endpoint_kernel(wave[n]))
  end

  endpoint_wave_inverse = ck_unit_series_inverse(
    endpoint_wave, order, identity_period, zero_period
  )

  slow_generator = typeof(identity_period)[zero_period for _ in 0:order]
  for n in 1:order
    slow_generator[n + 1] = ck_period_monomial(ck_endpoint_kernel(effective[n]), 1)
  end
  slow_propagator = ck_zero_constant_series_exponential(
    slow_generator, order, identity_period, zero_period
  )

  dressed_left = ck_truncated_series_product(
    endpoint_wave, slow_propagator, order, zero_period
  )
  amplitude = ck_truncated_series_product(
    dressed_left, endpoint_wave_inverse, order, zero_period
  )
  return CKPeriodAmplitudeResult(
    amplitude, endpoint_wave, endpoint_wave_inverse, slow_propagator
  )
end

function ck_endpoint_ordered_sideband_coefficient(
  state::CKEndpointKernel{H,T},
  output_channels::AbstractVector{Int},
  sidebands::AbstractVector{H};
  inverse_weight,
) where {H,T}
  length(output_channels) == length(sidebands) ||
    throw(ArgumentError("output channels and sidebands must have equal length"))

  result = state.zero_component
  for (key, value) in state.terms
    count(ck_is_jump, key.vertices) == length(sidebands) || continue

    prefix_mismatch = Vector{H}(undef, length(key.vertices) + 1)
    prefix_mismatch[1] = zero(H)
    sideband_index = 0
    valid = true
    for (vertex_index, vertex) in enumerate(key.vertices)
      mismatch = prefix_mismatch[vertex_index] + vertex.harmonic
      if ck_is_jump(vertex)
        sideband_index += 1
        if vertex.output_channel != output_channels[sideband_index]
          valid = false
          break
        end
        mismatch -= sidebands[sideband_index]
      end
      prefix_mismatch[vertex_index + 1] = mismatch
    end
    valid || continue

    weighted_value = value
    for interval in key.model_intervals
      interval_mismatch =
        prefix_mismatch[interval.stop + 1] - prefix_mismatch[interval.start + 1]
      if !iszero(interval_mismatch)
        valid = false
        break
      end
    end
    valid || continue

    for interval in key.resolvent_intervals
      interval_mismatch =
        prefix_mismatch[interval.stop + 1] - prefix_mismatch[interval.start + 1]
      if iszero(interval_mismatch)
        valid = false
        break
      end
      weighted_value = inverse_weight(interval_mismatch) * weighted_value
    end
    valid || continue
    result += weighted_value
  end
  return result
end

function ck_period_ordered_sideband_coefficients(
  polynomial::CKPeriodPolynomial{H,T},
  output_channels::AbstractVector{Int},
  sidebands::AbstractVector{H};
  inverse_weight,
) where {H,T}
  return T[
    ck_endpoint_ordered_sideband_coefficient(
      coefficient, output_channels, sidebands; inverse_weight=inverse_weight
    ) for coefficient in polynomial.coefficients
  ]
end

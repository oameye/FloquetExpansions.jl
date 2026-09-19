function ck_kernel_endpoint_product(
  left::CKPhysicalKernel{H,T}, right::CKPhysicalKernel{H,T}
) where {H,T}
  result = Dict{CKKernelKey{H},T}()
  for (left_key, left_value) in left.terms, (right_key, right_value) in right.terms
    right_length = length(right_key.vertices)
    vertices = CKBranchVertex{H}[right_key.vertices; left_key.vertices]
    model_cuts = ck_kernel_shifted_cuts(
      right_key.model_cuts, left_key.model_cuts, right_length
    )
    resolvent_cuts = ck_kernel_shifted_cuts(
      right_key.resolvent_cuts, left_key.resolvent_cuts, right_length
    )
    key = CKKernelKey(CKUnresolvedSector, vertices, model_cuts, resolvent_cuts)
    ck_kernel_accumulate!(result, key, left_value * right_value, left.zero_component)
  end
  return CKPhysicalKernel(result, left.zero_component)
end

function ck_kernel_zero_like(state::CKPhysicalKernel{H,T}) where {H,T}
  return CKPhysicalKernel(Dict{CKKernelKey{H},T}(), state.zero_component)
end

function ck_kernel_divide(state::CKPhysicalKernel{H,T}, divisor::Int) where {H,T}
  iszero(divisor) && throw(DivideError())
  result = Dict{CKKernelKey{H},T}()
  for (key, value) in state.terms
    ck_kernel_accumulate!(result, key, value / divisor, state.zero_component)
  end
  return CKPhysicalKernel(result, state.zero_component)
end

struct CKPeriodPolynomial{H,T}
  coefficients::Vector{CKPhysicalKernel{H,T}}
end

function ck_period_polynomial(coefficients::Vector{CKPhysicalKernel{H,T}}) where {H,T}
  isempty(coefficients) && throw(ArgumentError("period polynomial must contain a coefficient"))
  last_nonzero = findlast(coefficient -> !isempty(coefficient.terms), coefficients)
  length_to_keep = isnothing(last_nonzero) ? 1 : last_nonzero
  return CKPeriodPolynomial{H,T}(copy(coefficients[1:length_to_keep]))
end

function ck_period_zero(state::CKPhysicalKernel{H,T}) where {H,T}
  return ck_period_polynomial(CKPhysicalKernel{H,T}[ck_kernel_zero_like(state)])
end

function ck_period_constant(state::CKPhysicalKernel{H,T}) where {H,T}
  return ck_period_polynomial(CKPhysicalKernel{H,T}[state])
end

function ck_period_monomial(state::CKPhysicalKernel{H,T}, power::Int) where {H,T}
  power >= 0 || throw(ArgumentError("period power must be nonnegative"))
  coefficients = CKPhysicalKernel{H,T}[ck_kernel_zero_like(state) for _ in 0:power]
  coefficients[power + 1] = state
  return ck_period_polynomial(coefficients)
end

function Base.:(==)(left::CKPeriodPolynomial, right::CKPeriodPolynomial)
  return left.coefficients == right.coefficients
end

function Base.:+(left::CKPeriodPolynomial{H,T}, right::CKPeriodPolynomial{H,T}) where {H,T}
  zero_state = ck_kernel_zero_like(first(left.coefficients))
  length_result = max(length(left.coefficients), length(right.coefficients))
  coefficients = Vector{CKPhysicalKernel{H,T}}(undef, length_result)
  for index in 1:length_result
    left_value = index <= length(left.coefficients) ? left.coefficients[index] : zero_state
    right_value = index <= length(right.coefficients) ? right.coefficients[index] : zero_state
    coefficients[index] = left_value + right_value
  end
  return ck_period_polynomial(coefficients)
end

function Base.:-(left::CKPeriodPolynomial{H,T}, right::CKPeriodPolynomial{H,T}) where {H,T}
  zero_state = ck_kernel_zero_like(first(left.coefficients))
  length_result = max(length(left.coefficients), length(right.coefficients))
  coefficients = Vector{CKPhysicalKernel{H,T}}(undef, length_result)
  for index in 1:length_result
    left_value = index <= length(left.coefficients) ? left.coefficients[index] : zero_state
    right_value = index <= length(right.coefficients) ? right.coefficients[index] : zero_state
    coefficients[index] = left_value - right_value
  end
  return ck_period_polynomial(coefficients)
end

function Base.:-(value::CKPeriodPolynomial{H,T}) where {H,T}
  return ck_period_zero(first(value.coefficients)) - value
end

function Base.:*(left::CKPeriodPolynomial{H,T}, right::CKPeriodPolynomial{H,T}) where {H,T}
  zero_state = ck_kernel_zero_like(first(left.coefficients))
  maximum_power = length(left.coefficients) + length(right.coefficients) - 2
  coefficients = CKPhysicalKernel{H,T}[ck_kernel_zero_like(zero_state) for _ in 0:maximum_power]
  for left_index in eachindex(left.coefficients), right_index in eachindex(right.coefficients)
    output_index = left_index + right_index - 1
    coefficients[output_index] += ck_kernel_endpoint_product(
      left.coefficients[left_index], right.coefficients[right_index]
    )
  end
  return ck_period_polynomial(coefficients)
end

function Base.:/(value::CKPeriodPolynomial{H,T}, divisor::Int) where {H,T}
  return ck_period_polynomial([
    ck_kernel_divide(coefficient, divisor) for coefficient in value.coefficients
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
  series[1] == one_value || throw(ArgumentError("series must have unit order-zero coefficient"))

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

  zero_period = ck_period_zero(zero_state)
  identity_period = ck_period_constant(identity_state)

  endpoint_wave = typeof(identity_period)[zero_period for _ in 0:order]
  endpoint_wave[1] = identity_period
  for n in 1:(order - 1)
    endpoint_wave[n + 1] = ck_period_constant(wave[n])
  end

  endpoint_wave_inverse = ck_unit_series_inverse(
    endpoint_wave, order, identity_period, zero_period
  )

  slow_generator = typeof(identity_period)[zero_period for _ in 0:order]
  for n in 1:order
    slow_generator[n + 1] = ck_period_monomial(effective[n], 1)
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

function ck_period_ordered_sideband_coefficients(
  polynomial::CKPeriodPolynomial{H,T},
  output_channels::AbstractVector{Int},
  sidebands::AbstractVector{H};
  inverse_weight,
) where {H,T}
  return T[
    ck_kernel_ordered_sideband_coefficient(
      coefficient, output_channels, sidebands; inverse_weight=inverse_weight
    ) for coefficient in polynomial.coefficients
  ]
end

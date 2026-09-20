function ck_output_period_zero(zero_harmonic::H, zero_component::T) where {H,T}
  iszero(zero_harmonic) || throw(ArgumentError("zero harmonic must be zero"))
  kernel = CKOutputKernel(Dict{CKOutputKernelKey{H},T}(), zero_component)
  return CKOutputPeriodPolynomial(CKOutputKernel{H,T}[kernel])
end

function ck_output_period_add!(
  target::CKOutputPeriodPolynomial{H,T}, source::CKOutputPeriodPolynomial{H,T}
) where {H,T}
  while length(target.coefficients) < length(source.coefficients)
    push!(
      target.coefficients,
      CKOutputKernel(Dict{CKOutputKernelKey{H},T}(), target.coefficients[1].zero_component),
    )
  end
  for (period_index, source_kernel) in enumerate(source.coefficients)
    target_kernel = target.coefficients[period_index]
    for (key, value) in source_kernel.terms
      ck_output_accumulate!(
        target_kernel.terms, key, value, target_kernel.zero_component
      )
    end
  end
  return target
end

function ck_output_period_right_product(
  amplitude::CKOutputPeriodPolynomial{H,T},
  factor::CKOutputPairingPolynomial{T},
  zero_component::T,
) where {H,T}
  isempty(amplitude.coefficients) && throw(ArgumentError("amplitude polynomial must not be empty"))
  isempty(factor.terms) && return ck_output_period_zero(zero(H), zero_component)
  max_factor_power = maximum(keys(factor.terms))
  max_factor_power >= 0 || throw(ArgumentError("normalization contains a negative period power"))
  max_period_power = length(amplitude.coefficients) - 1 + max_factor_power
  coefficients = CKOutputKernel{H,T}[
    CKOutputKernel(Dict{CKOutputKernelKey{H},T}(), zero_component) for _ in 0:max_period_power
  ]

  for (amplitude_index, amplitude_kernel) in enumerate(amplitude.coefficients)
    amplitude_power = amplitude_index - 1
    for (key, amplitude_value) in amplitude_kernel.terms,
      (factor_power, factor_value) in factor.terms
      factor_power >= 0 || throw(ArgumentError("normalization contains a negative period power"))
      target = coefficients[amplitude_power + factor_power + 1]
      ck_output_accumulate!(
        target.terms, key, amplitude_value * factor_value, zero_component
      )
    end
  end
  return CKOutputPeriodPolynomial(coefficients)
end

function ck_output_right_normalize_series(
  amplitudes::Vector{CKOutputPeriodPolynomial{H,T}},
  normalization::CKOutputPairingSeries{T},
  order::Int,
  zero_component::T,
) where {H,T}
  order >= 0 || throw(ArgumentError("normalization order must be nonnegative"))
  result = CKOutputPeriodPolynomial{H,T}[
    ck_output_period_zero(zero(H), zero_component) for _ in 0:order
  ]

  for total_order in 0:order
    target = result[total_order + 1]
    for amplitude_order in 0:total_order
      amplitude_order + 1 <= length(amplitudes) || continue
      normalization_order = total_order - amplitude_order
      normalization_order + 1 <= length(normalization.coefficients) || continue
      term = ck_output_period_right_product(
        amplitudes[amplitude_order + 1],
        normalization.coefficients[normalization_order + 1],
        zero_component,
      )
      ck_output_period_add!(target, term)
    end
  end
  return result
end

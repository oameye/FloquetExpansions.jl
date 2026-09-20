struct CKOutputPairingSeries{H,T}
  coefficients::Vector{CKOutputPairingPolynomial{H,T}}
end

function ck_output_pairing_zero(::Type{H}, zero_component::T) where {H,T}
  return CKOutputPairingPolynomial(Dict{CKOutputPairingGrade{H},T}(), zero_component)
end

function ck_output_pairing_zero(zero_component::T) where {T}
  return ck_output_pairing_zero(Int, zero_component)
end

function ck_output_pairing_add!(
  target::CKOutputPairingPolynomial{H,T}, source::CKOutputPairingPolynomial{H,T}
) where {H,T}
  for (grade, value) in source.terms
    ck_output_pairing_accumulate!(target.terms, grade, value, target.zero_component)
  end
  return target
end

function ck_output_pairing_series(
  amplitudes::Vector{CKOutputPeriodPolynomial{H,S}},
  order::Int,
  imaginary::T,
  zero_component::R,
  operations::CKOutputPairingOperations,
) where {H<:Integer,S,T,R}
  order >= 0 || throw(ArgumentError("pairing order must be nonnegative"))
  result = [ck_output_pairing_zero(H, zero_component) for _ in 0:order]

  for total_order in 0:order
    for left_order in 0:total_order
      left_order + 1 <= length(amplitudes) || continue
      right_order = total_order - left_order
      right_order + 1 <= length(amplitudes) || continue
      contribution = ck_output_period_pairing(
        amplitudes[left_order + 1],
        amplitudes[right_order + 1],
        imaginary,
        zero_component,
        operations,
      )
      ck_output_pairing_add!(result[total_order + 1], contribution)
    end
  end
  return CKOutputPairingSeries(result)
end

function ck_output_channel_series(
  amplitudes::Vector{CKOutputPeriodPolynomial{H,S}},
  order::Int,
  imaginary::T,
  zero_component::R,
  channel_pair,
) where {H<:Integer,S,T,R}
  operations = CKOutputPairingOperations(identity, channel_pair, ck_output_channel_phase)
  return ck_output_pairing_series(amplitudes, order, imaginary, zero_component, operations)
end

function ck_output_metric_series(
  amplitudes::Vector{CKOutputPeriodPolynomial{H,S}},
  order::Int,
  imaginary::T,
  zero_component::R,
  metric_pair,
) where {H<:Integer,S,T,R}
  operations = CKOutputPairingOperations(conj, metric_pair, ck_output_metric_phase)
  return ck_output_pairing_series(amplitudes, order, imaginary, zero_component, operations)
end

function ck_output_pairing_series_coefficient(
  series::CKOutputPairingSeries, perturbative_order::Int
)
  perturbative_order >= 0 || throw(ArgumentError("perturbative order must be nonnegative"))
  perturbative_order + 1 <= length(series.coefficients) ||
    throw(BoundsError(series.coefficients, perturbative_order + 1))
  return series.coefficients[perturbative_order + 1]
end

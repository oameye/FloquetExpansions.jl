struct CKOutputPairingPolynomial{T}
  terms::Dict{Int,T}
  zero_component::T
end

function ck_output_pairing_accumulate!(
  terms::Dict{Int,T}, period_power::Int, value::T, zero_component::T
) where {T}
  updated = get(terms, period_power, zero_component) + value
  if ck_kernel_iszero(updated)
    haskey(terms, period_power) && delete!(terms, period_power)
  else
    terms[period_power] = updated
  end
  return terms
end

function ck_output_period_pairing(
  left::CKOutputPeriodPolynomial{H,S},
  right::CKOutputPeriodPolynomial{H,S},
  imaginary::T,
  zero_component::R,
  gram_weight,
  system_pair,
) where {H<:Integer,S,T,R}
  result = Dict{Int,R}()

  for (left_index, left_kernel) in enumerate(left.coefficients),
    (right_index, right_kernel) in enumerate(right.coefficients)
    left_period_power = left_index - 1
    right_period_power = right_index - 1

    for (left_key, left_value) in left_kernel.terms,
      (right_key, right_value) in right_kernel.terms
      left_key.output_channels == right_key.output_channels || continue
      output_number = length(left_key.output_channels)
      inverse_fourier_shift = left_period_power + right_period_power - 2 * output_number
      overlap = ck_output_time_overlap(left_key, right_key, imaginary)
      paired_system = system_pair(left_value, right_value)

      for (overlap_period_power, overlap_value) in overlap.terms
        total_period_power = inverse_fourier_shift + overlap_period_power
        value = gram_weight(overlap_value) * paired_system
        ck_output_pairing_accumulate!(result, total_period_power, value, zero_component)
      end
    end
  end
  return CKOutputPairingPolynomial(result, zero_component)
end

function ck_output_channel_pairing(
  left::CKOutputPeriodPolynomial{H,S},
  right::CKOutputPeriodPolynomial{H,S},
  imaginary::T,
  zero_component::R,
  channel_pair,
) where {H<:Integer,S,T,R}
  return ck_output_period_pairing(
    left, right, imaginary, zero_component, identity, channel_pair
  )
end

function ck_output_metric_pairing(
  left::CKOutputPeriodPolynomial{H,S},
  right::CKOutputPeriodPolynomial{H,S},
  imaginary::T,
  zero_component::R,
  metric_pair,
) where {H<:Integer,S,T,R}
  return ck_output_period_pairing(left, right, imaginary, zero_component, conj, metric_pair)
end

function ck_output_pairing_coefficient(
  polynomial::CKOutputPairingPolynomial{T}, period_power::Int
) where {T}
  return get(polynomial.terms, period_power, polynomial.zero_component)
end

function ck_output_pairing_has_negative_power(polynomial::CKOutputPairingPolynomial)
  return any(period_power -> period_power < 0, keys(polynomial.terms))
end

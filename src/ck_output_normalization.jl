function ck_output_pairing_identity(identity_component::T) where {T}
  return CKOutputPairingPolynomial(Dict(0 => identity_component), zero(identity_component))
end

function ck_output_pairing_product(
  left::CKOutputPairingPolynomial{T},
  right::CKOutputPairingPolynomial{T},
  zero_component::T,
) where {T}
  result = Dict{Int,T}()
  for (left_power, left_value) in left.terms, (right_power, right_value) in right.terms
    ck_output_pairing_accumulate!(
      result, left_power + right_power, left_value * right_value, zero_component
    )
  end
  return CKOutputPairingPolynomial(result, zero_component)
end

function ck_output_pairing_scale(
  scale, polynomial::CKOutputPairingPolynomial{T}, zero_component::T
) where {T}
  result = Dict{Int,T}()
  for (period_power, value) in polynomial.terms
    ck_output_pairing_accumulate!(result, period_power, scale * value, zero_component)
  end
  return CKOutputPairingPolynomial(result, zero_component)
end

function ck_output_pairing_series_zero(order::Int, zero_component::T) where {T}
  order >= 0 || throw(ArgumentError("series order must be nonnegative"))
  return CKOutputPairingSeries([
    ck_output_pairing_zero(zero_component) for _ in 0:order
  ])
end

function ck_output_pairing_series_product(
  left::CKOutputPairingSeries{T},
  right::CKOutputPairingSeries{T},
  order::Int,
  zero_component::T,
) where {T}
  order >= 0 || throw(ArgumentError("series order must be nonnegative"))
  result = ck_output_pairing_series_zero(order, zero_component)
  for total_order in 0:order
    target = result.coefficients[total_order + 1]
    for left_order in 0:total_order
      left_order + 1 <= length(left.coefficients) || continue
      right_order = total_order - left_order
      right_order + 1 <= length(right.coefficients) || continue
      term = ck_output_pairing_product(
        left.coefficients[left_order + 1], right.coefficients[right_order + 1], zero_component
      )
      ck_output_pairing_add!(target, term)
    end
  end
  return result
end

function ck_output_metric_inverse_sqrt_series(
  metric::CKOutputPairingSeries{T}, order::Int, identity_component::T
) where {T}
  order >= 0 || throw(ArgumentError("normalization order must be nonnegative"))
  isempty(metric.coefficients) && throw(ArgumentError("metric series must contain a leading coefficient"))
  zero_component = zero(identity_component)
  leading = metric.coefficients[1]
  leading.terms == Dict(0 => identity_component) ||
    throw(ArgumentError("metric series must have unit leading coefficient"))

  coefficients = [ck_output_pairing_zero(zero_component) for _ in 0:order]
  coefficients[1] = ck_output_pairing_identity(identity_component)
  normalization = CKOutputPairingSeries(coefficients)

  for total_order in 1:order
    residual = ck_output_pairing_zero(zero_component)
    for left_order in 0:total_order
      left_order + 1 <= length(normalization.coefficients) || continue
      for metric_order in 0:(total_order - left_order)
        metric_order + 1 <= length(metric.coefficients) || continue
        right_order = total_order - left_order - metric_order
        right_order + 1 <= length(normalization.coefficients) || continue
        left_metric = ck_output_pairing_product(
          normalization.coefficients[left_order + 1],
          metric.coefficients[metric_order + 1],
          zero_component,
        )
        term = ck_output_pairing_product(
          left_metric, normalization.coefficients[right_order + 1], zero_component
        )
        ck_output_pairing_add!(residual, term)
      end
    end
    normalization.coefficients[total_order + 1] = ck_output_pairing_scale(
      -1 // 2, residual, zero_component
    )
  end
  return normalization
end

function ck_output_normalized_metric_series(
  metric::CKOutputPairingSeries{T},
  normalization::CKOutputPairingSeries{T},
  order::Int,
  zero_component::T,
) where {T}
  left = ck_output_pairing_series_product(normalization, metric, order, zero_component)
  return ck_output_pairing_series_product(left, normalization, order, zero_component)
end

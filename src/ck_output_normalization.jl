function ck_output_pairing_identity(zero_harmonic::H, identity_component::T) where {H,T}
  grade = CKOutputPairingGrade(zero_harmonic, 0)
  return CKOutputPairingPolynomial(Dict(grade => identity_component), zero(identity_component))
end

function ck_output_pairing_identity(identity_component::T) where {T}
  return ck_output_pairing_identity(0, identity_component)
end

function ck_output_pairing_product(
  left::CKOutputPairingPolynomial{H,T},
  right::CKOutputPairingPolynomial{H,T},
  zero_component::T,
) where {H,T}
  result = Dict{CKOutputPairingGrade{H},T}()
  for (left_grade, left_value) in left.terms, (right_grade, right_value) in right.terms
    grade = CKOutputPairingGrade(
      left_grade.phase_harmonic + right_grade.phase_harmonic,
      left_grade.period_power + right_grade.period_power,
    )
    ck_output_pairing_accumulate!(result, grade, left_value * right_value, zero_component)
  end
  return CKOutputPairingPolynomial(result, zero_component)
end

function ck_output_pairing_scale(
  scale, polynomial::CKOutputPairingPolynomial{H,T}, zero_component::T
) where {H,T}
  result = Dict{CKOutputPairingGrade{H},T}()
  for (grade, value) in polynomial.terms
    ck_output_pairing_accumulate!(result, grade, scale * value, zero_component)
  end
  return CKOutputPairingPolynomial(result, zero_component)
end

function ck_output_pairing_series_zero(
  order::Int, zero_harmonic::H, zero_component::T
) where {H,T}
  order >= 0 || throw(ArgumentError("series order must be nonnegative"))
  return CKOutputPairingSeries([
    ck_output_pairing_zero(zero_harmonic, zero_component) for _ in 0:order
  ])
end

function ck_output_pairing_series_zero(order::Int, zero_component::T) where {T}
  return ck_output_pairing_series_zero(order, 0, zero_component)
end

function ck_output_pairing_series_product(
  left::CKOutputPairingSeries{H,T},
  right::CKOutputPairingSeries{H,T},
  order::Int,
  zero_component::T,
) where {H,T}
  order >= 0 || throw(ArgumentError("series order must be nonnegative"))
  result = ck_output_pairing_series_zero(order, zero(H), zero_component)
  for total_order in 0:order
    target = result.coefficients[total_order + 1]
    for left_order in 0:total_order
      left_order + 1 <= length(left.coefficients) || continue
      right_order = total_order - left_order
      right_order + 1 <= length(right.coefficients) || continue
      term = ck_output_pairing_product(
        left.coefficients[left_order + 1],
        right.coefficients[right_order + 1],
        zero_component,
      )
      ck_output_pairing_add!(target, term)
    end
  end
  return result
end

function ck_output_pairing_is_identity(
  polynomial::CKOutputPairingPolynomial,
  grade::CKOutputPairingGrade,
  identity_component,
)
  length(polynomial.terms) == 1 || return false
  for (key, value) in polynomial.terms
    return isequal(key, grade) && isequal(value, identity_component)
  end
  return false
end

function ck_output_metric_inverse_sqrt_residual(
  metric::CKOutputPairingSeries{H,T},
  normalization::CKOutputPairingSeries{H,T},
  total_order::Int,
  zero_component::T,
) where {H,T}
  residual = ck_output_pairing_zero(zero(H), zero_component)
  maximum_metric_order = length(metric.coefficients) - 1
  for left_order in 0:total_order
    upper_metric_order = min(total_order - left_order, maximum_metric_order)
    for metric_order in 0:upper_metric_order
      right_order = total_order - left_order - metric_order
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
  return residual
end

function ck_output_metric_inverse_sqrt_series(
  metric::CKOutputPairingSeries{H,T}, order::Int, identity_component::T
) where {H,T}
  order >= 0 || throw(ArgumentError("normalization order must be nonnegative"))
  isempty(metric.coefficients) &&
    throw(ArgumentError("metric series must contain a leading coefficient"))
  zero_component = zero(identity_component)
  identity_grade = CKOutputPairingGrade(zero(H), 0)
  ck_output_pairing_is_identity(metric.coefficients[1], identity_grade, identity_component) ||
    throw(ArgumentError("metric series must have unit leading coefficient"))

  coefficients = [ck_output_pairing_zero(zero(H), zero_component) for _ in 0:order]
  coefficients[1] = ck_output_pairing_identity(zero(H), identity_component)
  normalization = CKOutputPairingSeries(coefficients)

  for total_order in 1:order
    residual = ck_output_metric_inverse_sqrt_residual(
      metric, normalization, total_order, zero_component
    )
    normalization.coefficients[total_order + 1] = ck_output_pairing_scale(
      -1 // 2, residual, zero_component
    )
  end
  return normalization
end

function ck_output_normalized_metric_series(
  metric::CKOutputPairingSeries{H,T},
  normalization::CKOutputPairingSeries{H,T},
  order::Int,
  zero_component::T,
) where {H,T}
  left = ck_output_pairing_series_product(normalization, metric, order, zero_component)
  return ck_output_pairing_series_product(left, normalization, order, zero_component)
end

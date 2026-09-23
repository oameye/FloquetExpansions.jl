struct CKOutputPairingGrade{H}
  phase_harmonic::H
  period_power::Int
end

function Base.isequal(left::CKOutputPairingGrade, right::CKOutputPairingGrade)
  return isequal(left.phase_harmonic, right.phase_harmonic) &&
         left.period_power == right.period_power
end

Base.:(==)(left::CKOutputPairingGrade, right::CKOutputPairingGrade) = isequal(left, right)

function Base.hash(grade::CKOutputPairingGrade, seed::UInt)
  return hash(grade.period_power, hash(grade.phase_harmonic, seed))
end

struct CKOutputPairingPolynomial{H,T}
  terms::Dict{CKOutputPairingGrade{H},T}
  zero_component::T
end

function CKOutputPairingPolynomial(terms::Dict{Int,T}, zero_component::T) where {T}
  graded_terms = Dict{CKOutputPairingGrade{Int},T}()
  for (period_power, value) in terms
    graded_terms[CKOutputPairingGrade(0, period_power)] = value
  end
  return CKOutputPairingPolynomial(graded_terms, zero_component)
end

struct CKOutputPairingOperations{G,S,P}
  gram_weight::G
  system_pair::S
  phase_pair::P
end

function Base.isequal(left::CKOutputPairingPolynomial, right::CKOutputPairingPolynomial)
  return isequal(left.zero_component, right.zero_component) &&
         isequal(left.terms, right.terms)
end

function Base.:(==)(left::CKOutputPairingPolynomial, right::CKOutputPairingPolynomial)
  return isequal(left, right)
end

function ck_output_pairing_accumulate!(
  terms::Dict{CKOutputPairingGrade{H},T},
  grade::CKOutputPairingGrade{H},
  value::T,
  zero_component::T,
) where {H,T}
  updated = get(terms, grade, zero_component) + value
  if ck_kernel_iszero(updated)
    haskey(terms, grade) && delete!(terms, grade)
  else
    terms[grade] = updated
  end
  return terms
end

function ck_output_pairing_accumulate!(
  terms::Dict{CKOutputPairingGrade{H},T},
  phase_harmonic::H,
  period_power::Int,
  value::T,
  zero_component::T,
) where {H,T}
  return ck_output_pairing_accumulate!(
    terms, CKOutputPairingGrade(phase_harmonic, period_power), value, zero_component
  )
end

ck_output_channel_phase(left, right) = left - right
ck_output_metric_phase(left, right) = right - left

function ck_output_period_pairing(
  left::CKOutputPeriodPolynomial{H,S},
  right::CKOutputPeriodPolynomial{H,S},
  imaginary::T,
  zero_component::R,
  operations::CKOutputPairingOperations,
) where {H<:Integer,S,T,R}
  result = Dict{CKOutputPairingGrade{H},R}()

  for (left_index, left_kernel) in enumerate(left.coefficients),
    (right_index, right_kernel) in enumerate(right.coefficients)

    left_period_power = left_index - 1
    right_period_power = right_index - 1

    for (left_key, left_value) in left_kernel.terms,
      (right_key, right_value) in right_kernel.terms

      isequal(left_key.output_channels, right_key.output_channels) || continue
      output_number = length(left_key.output_channels)
      inverse_fourier_shift = left_period_power + right_period_power - 2 * output_number
      overlap = ck_output_time_overlap(left_key, right_key, imaginary)
      paired_system = operations.system_pair(left_value, right_value)
      phase_harmonic = operations.phase_pair(
        left_key.phase_harmonic, right_key.phase_harmonic
      )

      for (overlap_period_power, overlap_value) in overlap.terms
        total_period_power = inverse_fourier_shift + overlap_period_power
        value = operations.gram_weight(overlap_value) * paired_system
        ck_output_pairing_accumulate!(
          result, phase_harmonic, total_period_power, value, zero_component
        )
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
  operations = CKOutputPairingOperations(identity, channel_pair, ck_output_channel_phase)
  return ck_output_period_pairing(left, right, imaginary, zero_component, operations)
end

function ck_output_metric_pairing(
  left::CKOutputPeriodPolynomial{H,S},
  right::CKOutputPeriodPolynomial{H,S},
  imaginary::T,
  zero_component::R,
  metric_pair,
) where {H<:Integer,S,T,R}
  operations = CKOutputPairingOperations(conj, metric_pair, ck_output_metric_phase)
  return ck_output_period_pairing(left, right, imaginary, zero_component, operations)
end

function ck_output_pairing_coefficient(
  polynomial::CKOutputPairingPolynomial{H,T}, phase_harmonic::H, period_power::Int
) where {H,T}
  grade = CKOutputPairingGrade(phase_harmonic, period_power)
  return get(polynomial.terms, grade, polynomial.zero_component)
end

function ck_output_pairing_coefficient(
  polynomial::CKOutputPairingPolynomial{H,T}, period_power::Int
) where {H,T}
  result = polynomial.zero_component
  for (grade, value) in polynomial.terms
    grade.period_power == period_power || continue
    result += value
  end
  return result
end

function ck_output_pairing_period_terms(
  polynomial::CKOutputPairingPolynomial{H,T}
) where {H,T}
  result = Dict{Int,T}()
  for (grade, value) in polynomial.terms
    updated = get(result, grade.period_power, polynomial.zero_component) + value
    if ck_kernel_iszero(updated)
      haskey(result, grade.period_power) && delete!(result, grade.period_power)
    else
      result[grade.period_power] = updated
    end
  end
  return result
end

function ck_output_pairing_phase_support(polynomial::CKOutputPairingPolynomial{H}) where {H}
  phases = H[]
  for grade in keys(polynomial.terms)
    grade.phase_harmonic in phases || push!(phases, grade.phase_harmonic)
  end
  sort!(phases)
  return phases
end

function ck_output_pairing_has_negative_power(polynomial::CKOutputPairingPolynomial)
  return any(grade -> grade.period_power < 0, keys(polynomial.terms))
end

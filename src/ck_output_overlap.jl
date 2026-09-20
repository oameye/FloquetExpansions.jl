struct CKOutputTimeMonomial{H,T}
  frequencies::Vector{H}
  degrees::Vector{Int}
  period_power::Int
  coefficient::T
end

struct CKOutputOverlapPolynomial{T}
  terms::Dict{Int,T}
  zero_coefficient::T
end

function ck_output_overlap_accumulate!(
  terms::Dict{Int,T}, period_power::Int, value::T, zero_coefficient::T
) where {T}
  updated = get(terms, period_power, zero_coefficient) + value
  if iszero(updated)
    haskey(terms, period_power) && delete!(terms, period_power)
  else
    terms[period_power] = updated
  end
  return terms
end

function ck_output_time_monomial_product(
  left::CKOutputTimeMonomial{H,T}, right::CKOutputTimeMonomial{H,T}
) where {H,T}
  length(left.frequencies) == length(right.frequencies) ||
    throw(ArgumentError("time monomials must have equal output number"))
  frequencies = H[
    left.frequencies[index] + right.frequencies[index] for
    index in eachindex(left.frequencies)
  ]
  degrees = Int[
    left.degrees[index] + right.degrees[index] for index in eachindex(left.degrees)
  ]
  return CKOutputTimeMonomial(
    frequencies,
    degrees,
    left.period_power + right.period_power,
    left.coefficient * right.coefficient,
  )
end

function ck_output_coordinate_frequencies(
  block::CKOutputBlock, output_stop::Int, harmonic::H, output_number::Int
) where {H<:Integer}
  frequencies = fill(zero(H), output_number)
  frequencies[output_stop] += harmonic
  if output_stop < block.output_stop
    frequencies[output_stop + 1] -= harmonic
  end
  return frequencies
end

function ck_output_gap_polynomial_terms(
  block::CKOutputBlock,
  output_stop::Int,
  output_number::Int,
  harmonic::H,
  period_power::Int,
  coordinate_power::Int,
  coefficient::T,
) where {H<:Integer,T}
  frequencies = ck_output_coordinate_frequencies(
    block, output_stop, harmonic, output_number
  )
  result = CKOutputTimeMonomial{H,T}[]

  if output_stop == block.output_stop
    degrees = zeros(Int, output_number)
    degrees[output_stop] = coordinate_power
    push!(result, CKOutputTimeMonomial(frequencies, degrees, period_power, coefficient))
    return result
  end

  # For tau_j < tau_{j+1}, the periodic cumulative-sideband coordinate is
  # x_j = L + tau_j - tau_{j+1}.  Expand its finite polynomial exactly.
  for difference_power in 0:coordinate_power
    period_increment = coordinate_power - difference_power
    choose_period = binomial(coordinate_power, difference_power)
    for left_power in 0:difference_power
      right_power = difference_power - left_power
      sign = isodd(right_power) ? -1 : 1
      value = coefficient * choose_period * binomial(difference_power, left_power) * sign
      iszero(value) && continue
      degrees = zeros(Int, output_number)
      degrees[output_stop] = left_power
      degrees[output_stop + 1] = right_power
      push!(
        result,
        CKOutputTimeMonomial(
          copy(frequencies), degrees, period_power + period_increment, value
        ),
      )
    end
  end
  return result
end

function ck_output_coordinate_time_terms(
  coordinate::CKOutputTimeCoordinate{H,T}, block::CKOutputBlock, output_number::Int
) where {H<:Integer,T}
  result = CKOutputTimeMonomial{H,T}[]

  for phase in coordinate.kernel.phases
    frequencies = ck_output_coordinate_frequencies(
      block, coordinate.output_stop, phase.harmonic, output_number
    )
    push!(
      result,
      CKOutputTimeMonomial(frequencies, zeros(Int, output_number), 0, phase.coefficient),
    )
  end

  for component in coordinate.kernel.homological
    polynomial = ck_homological_kernel(component.order, one(component.coefficient))
    for ((period_power, coordinate_power), value) in polynomial.terms
      coefficient = component.coefficient * value
      append!(
        result,
        ck_output_gap_polynomial_terms(
          block,
          coordinate.output_stop,
          output_number,
          component.harmonic,
          period_power,
          coordinate_power,
          coefficient,
        ),
      )
    end
  end
  return result
end

function ck_output_time_monomials(
  key::CKOutputKernelKey{H}, realization::CKOutputTimeRealization{H,T}
) where {H<:Integer,T}
  realization.valid || return CKOutputTimeMonomial{H,T}[]
  output_number = length(key.output_channels)
  terms = CKOutputTimeMonomial{H,T}[CKOutputTimeMonomial(
    fill(zero(H), output_number), zeros(Int, output_number), 0, realization.scalar
  ),]

  for coordinate in realization.coordinates
    block = realization.blocks[coordinate.block]
    coordinate_terms = ck_output_coordinate_time_terms(coordinate, block, output_number)
    next_terms = CKOutputTimeMonomial{H,T}[]
    for term in terms, coordinate_term in coordinate_terms
      push!(next_terms, ck_output_time_monomial_product(term, coordinate_term))
    end
    terms = next_terms
  end
  return terms
end

function ck_output_simplex_state_accumulate!(
  state::Dict{H,Vector{T}}, frequency::H, degree::Int, value::T, zero_coefficient::T
) where {H<:Integer,T}
  coefficients = get!(state, frequency) do
    return T[]
  end
  while length(coefficients) <= degree
    push!(coefficients, zero_coefficient)
  end
  coefficients[degree + 1] += value
  return state
end

function ck_output_simplex_state_scale_add!(
  target::Dict{H,Vector{T}}, source::Dict{H,Vector{T}}, scale::T, zero_coefficient::T
) where {H<:Integer,T}
  for (frequency, coefficients) in source
    for (degree_index, coefficient) in enumerate(coefficients)
      ck_output_simplex_state_accumulate!(
        target, frequency, degree_index - 1, scale * coefficient, zero_coefficient
      )
    end
  end
  return target
end

function ck_output_simplex_monomial_primitive(
  frequency::H, degree::Int, imaginary::T
) where {H<:Integer,T}
  degree >= 0 || throw(ArgumentError("simplex monomial degree must be nonnegative"))
  zero_coefficient = zero(imaginary)
  result = Dict{H,Vector{T}}()

  if iszero(frequency)
    ck_output_simplex_state_accumulate!(
      result, zero(H), degree + 1, one(imaginary) / (degree + 1), zero_coefficient
    )
    return result
  end

  inverse_frequency = one(imaginary) / (imaginary * frequency)
  if iszero(degree)
    ck_output_simplex_state_accumulate!(
      result, zero(H), 0, inverse_frequency, zero_coefficient
    )
    ck_output_simplex_state_accumulate!(
      result, frequency, 0, -inverse_frequency, zero_coefficient
    )
    return result
  end

  lower = ck_output_simplex_monomial_primitive(frequency, degree - 1, imaginary)
  ck_output_simplex_state_scale_add!(
    result, lower, degree * inverse_frequency, zero_coefficient
  )
  ck_output_simplex_state_accumulate!(
    result, frequency, degree, -inverse_frequency, zero_coefficient
  )
  return result
end

function ck_output_simplex_shift_state(
  state::Dict{H,Vector{T}}, frequency::H, degree::Int, zero_coefficient::T
) where {H<:Integer,T}
  result = Dict{H,Vector{T}}()
  for (state_frequency, coefficients) in state
    for (degree_index, coefficient) in enumerate(coefficients)
      ck_output_simplex_state_accumulate!(
        result,
        state_frequency + frequency,
        degree_index - 1 + degree,
        coefficient,
        zero_coefficient,
      )
    end
  end
  return result
end

function ck_output_simplex_integrate_state(
  state::Dict{H,Vector{T}}, imaginary::T
) where {H<:Integer,T}
  result = Dict{H,Vector{T}}()
  zero_coefficient = zero(imaginary)
  for (frequency, coefficients) in state
    for (degree_index, coefficient) in enumerate(coefficients)
      primitive = ck_output_simplex_monomial_primitive(
        frequency, degree_index - 1, imaginary
      )
      ck_output_simplex_state_scale_add!(result, primitive, coefficient, zero_coefficient)
    end
  end
  return result
end

function ck_output_simplex_monomial_integral(
  frequencies::AbstractVector{H}, degrees::AbstractVector{Int}, imaginary::T
) where {H<:Integer,T}
  length(frequencies) == length(degrees) ||
    throw(ArgumentError("simplex frequencies and degrees must have equal length"))
  zero_coefficient = zero(imaginary)
  state = Dict{H,Vector{T}}(zero(H) => T[one(imaginary)])

  for index in eachindex(frequencies)
    state = ck_output_simplex_shift_state(
      state, frequencies[index], degrees[index], zero_coefficient
    )
    state = ck_output_simplex_integrate_state(state, imaginary)
  end

  terms = Dict{Int,T}()
  for coefficients in values(state)
    for (degree_index, coefficient) in enumerate(coefficients)
      ck_output_overlap_accumulate!(terms, degree_index - 1, coefficient, zero_coefficient)
    end
  end
  return CKOutputOverlapPolynomial(terms, zero_coefficient)
end

function ck_output_time_overlap(
  left_key::CKOutputKernelKey{H}, right_key::CKOutputKernelKey{H}, imaginary::T
) where {H<:Integer,T}
  zero_coefficient = zero(imaginary)
  left_key.output_channels == right_key.output_channels ||
    return CKOutputOverlapPolynomial(Dict{Int,T}(), zero_coefficient)

  left_realization = ck_output_time_realization(left_key, imaginary)
  right_realization = ck_output_time_realization(right_key, imaginary)
  left_realization.valid && right_realization.valid ||
    return CKOutputOverlapPolynomial(Dict{Int,T}(), zero_coefficient)

  left_terms = ck_output_time_monomials(left_key, left_realization)
  right_terms = ck_output_time_monomials(right_key, right_realization)
  overlap_terms = Dict{Int,T}()

  for left_term in left_terms, right_term in right_terms
    frequencies = H[
      left_term.frequencies[index] - right_term.frequencies[index] for
      index in eachindex(left_term.frequencies)
    ]
    degrees = Int[
      left_term.degrees[index] + right_term.degrees[index] for
      index in eachindex(left_term.degrees)
    ]
    simplex = ck_output_simplex_monomial_integral(frequencies, degrees, imaginary)
    scale = left_term.coefficient * conj(right_term.coefficient)
    period_shift = left_term.period_power + right_term.period_power
    for (period_power, coefficient) in simplex.terms
      ck_output_overlap_accumulate!(
        overlap_terms, period_shift + period_power, scale * coefficient, zero_coefficient
      )
    end
  end
  return CKOutputOverlapPolynomial(overlap_terms, zero_coefficient)
end

function ck_output_overlap_coefficient(
  polynomial::CKOutputOverlapPolynomial{T}, period_power::Int
) where {T}
  period_power >= 0 || throw(ArgumentError("period power must be nonnegative"))
  return get(polynomial.terms, period_power, polynomial.zero_coefficient)
end

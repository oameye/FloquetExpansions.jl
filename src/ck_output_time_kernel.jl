struct CKHomologicalPolynomial{T}
  terms::Dict{Tuple{Int,Int},T}
  zero_coefficient::T
end

function ck_time_polynomial_accumulate!(
  terms::Dict{Tuple{Int,Int},T},
  period_power::Int,
  coordinate_power::Int,
  value::T,
  zero_coefficient::T,
) where {T}
  key = (period_power, coordinate_power)
  updated = get(terms, key, zero_coefficient) + value
  if iszero(updated)
    haskey(terms, key) && delete!(terms, key)
  else
    terms[key] = updated
  end
  return terms
end

function ck_homological_polynomial(
  terms::Dict{Tuple{Int,Int},T}, zero_coefficient::T
) where {T}
  result = Dict{Tuple{Int,Int},T}()
  for ((period_power, coordinate_power), value) in terms
    period_power >= 0 || throw(ArgumentError("period power must be nonnegative"))
    coordinate_power >= 0 || throw(ArgumentError("coordinate power must be nonnegative"))
    iszero(value) ||
      ck_time_polynomial_accumulate!(
        result, period_power, coordinate_power, value, zero_coefficient
      )
  end
  return CKHomologicalPolynomial(result, zero_coefficient)
end

function ck_homological_average(polynomial::CKHomologicalPolynomial{T}) where {T}
  result = Dict{Int,T}()
  for ((period_power, coordinate_power), value) in polynomial.terms
    total_power = period_power + coordinate_power
    updated = get(result, total_power, polynomial.zero_coefficient) +
              value / (coordinate_power + 1)
    if iszero(updated)
      haskey(result, total_power) && delete!(result, total_power)
    else
      result[total_power] = updated
    end
  end
  return result
end

function ck_homological_next(polynomial::CKHomologicalPolynomial{T}) where {T}
  terms = Dict{Tuple{Int,Int},T}()
  for ((period_power, coordinate_power), value) in polynomial.terms
    ck_time_polynomial_accumulate!(
      terms,
      period_power,
      coordinate_power + 1,
      -value / (coordinate_power + 1),
      polynomial.zero_coefficient,
    )
  end

  primitive = ck_homological_polynomial(terms, polynomial.zero_coefficient)
  for (period_power, average) in ck_homological_average(primitive)
    ck_time_polynomial_accumulate!(
      terms, period_power, 0, -average, polynomial.zero_coefficient
    )
  end
  return ck_homological_polynomial(terms, polynomial.zero_coefficient)
end

function ck_homological_kernel(order::Int, one_coefficient::T) where {T}
  order >= 1 || throw(ArgumentError("homological order must be positive"))
  zero_coefficient = zero(one_coefficient)
  terms = Dict{Tuple{Int,Int},T}(
    (0, 1) => one_coefficient,
    (1, 0) => -one_coefficient / 2,
  )
  result = ck_homological_polynomial(terms, zero_coefficient)
  for _ in 2:order
    result = ck_homological_next(result)
  end
  return result
end

function ck_homological_derivative(polynomial::CKHomologicalPolynomial{T}) where {T}
  terms = Dict{Tuple{Int,Int},T}()
  for ((period_power, coordinate_power), value) in polynomial.terms
    iszero(coordinate_power) && continue
    ck_time_polynomial_accumulate!(
      terms,
      period_power,
      coordinate_power - 1,
      coordinate_power * value,
      polynomial.zero_coefficient,
    )
  end
  return ck_homological_polynomial(terms, polynomial.zero_coefficient)
end

function ck_homological_scale(polynomial::CKHomologicalPolynomial{T}, scale::T) where {T}
  terms = Dict{Tuple{Int,Int},T}()
  for ((period_power, coordinate_power), value) in polynomial.terms
    ck_time_polynomial_accumulate!(
      terms,
      period_power,
      coordinate_power,
      scale * value,
      polynomial.zero_coefficient,
    )
  end
  return ck_homological_polynomial(terms, polynomial.zero_coefficient)
end

struct CKHomologicalComponent{H,T}
  harmonic::H
  order::Int
  coefficient::T
end

struct CKPhaseComponent{H,T}
  harmonic::H
  coefficient::T
end

struct CKResolventTimeKernel{H,T}
  homological::Vector{CKHomologicalComponent{H,T}}
  phases::Vector{CKPhaseComponent{H,T}}
  zero_coefficient::T
end

function ck_resolvent_pole_multiplicities(poles::AbstractVector{H}) where {H<:Integer}
  multiplicities = Dict{H,Int}()
  for pole in poles
    multiplicities[pole] = get(multiplicities, pole, 0) + 1
  end
  return multiplicities
end

function ck_resolvent_taylor_factor(
  difference::H, multiplicity::Int, max_degree::Int, one_coefficient::T
) where {H<:Integer,T}
  iszero(difference) && throw(ArgumentError("partial-fraction centers must be distinct"))
  result = Vector{T}(undef, max_degree + 1)
  for degree in 0:max_degree
    coefficient =
      one_coefficient * binomial(multiplicity + degree - 1, degree) /
      difference^(multiplicity + degree)
    result[degree + 1] = isodd(degree) ? -coefficient : coefficient
  end
  return result
end

function ck_resolvent_series_product(
  left::Vector{T}, right::Vector{T}, zero_coefficient::T
) where {T}
  max_degree = length(left) - 1
  length(right) == length(left) || throw(ArgumentError("series lengths must agree"))
  result = fill(zero_coefficient, max_degree + 1)
  for left_degree in 0:max_degree, right_degree in 0:(max_degree - left_degree)
    result[left_degree + right_degree + 1] +=
      left[left_degree + 1] * right[right_degree + 1]
  end
  return result
end

function ck_resolvent_time_kernel(
  poles::AbstractVector{H}, imaginary::T
) where {H<:Integer,T}
  isempty(poles) && throw(ArgumentError("at least one resolvent pole is required"))
  zero_coefficient = zero(imaginary)
  one_coefficient = one(imaginary)
  multiplicities = ck_resolvent_pole_multiplicities(poles)
  centers = sort!(collect(keys(multiplicities)))
  total_order = length(poles)
  homological = CKHomologicalComponent{H,T}[]

  for center in centers
    multiplicity = multiplicities[center]
    max_degree = multiplicity - 1
    series = fill(zero_coefficient, max_degree + 1)
    series[1] = one_coefficient

    for other_center in centers
      other_center == center && continue
      factor = ck_resolvent_taylor_factor(
        other_center - center,
        multiplicities[other_center],
        max_degree,
        one_coefficient,
      )
      series = ck_resolvent_series_product(series, factor, zero_coefficient)
    end

    for order in 1:multiplicity
      taylor_degree = multiplicity - order
      coefficient = imaginary^(total_order - order) * series[taylor_degree + 1]
      iszero(coefficient) ||
        push!(homological, CKHomologicalComponent(center, order, coefficient))
    end
  end

  phases = CKPhaseComponent{H,T}[]
  for center in centers
    residual = zero_coefficient
    for component in homological
      component.harmonic == center && continue
      mismatch = component.harmonic - center
      residual += component.coefficient * (imaginary / mismatch)^component.order
    end
    iszero(residual) || push!(phases, CKPhaseComponent(center, -residual))
  end

  return CKResolventTimeKernel(homological, phases, zero_coefficient)
end

function ck_resolvent_time_fourier_coefficient(
  kernel::CKResolventTimeKernel{H,T}, sideband::H, imaginary::T
) where {H<:Integer,T}
  result = kernel.zero_coefficient
  for phase in kernel.phases
    phase.harmonic == sideband && (result += phase.coefficient)
  end
  for component in kernel.homological
    mismatch = component.harmonic - sideband
    iszero(mismatch) && continue
    result += component.coefficient * (imaginary / mismatch)^component.order
  end
  return result
end

function ck_resolvent_product_fourier_coefficient(
  poles::AbstractVector{H}, sideband::H, imaginary::T
) where {H<:Integer,T}
  result = one(imaginary)
  for pole in poles
    mismatch = pole - sideband
    iszero(mismatch) && return zero(imaginary)
    result *= imaginary / mismatch
  end
  return result
end

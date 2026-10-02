function canonical_liouvillian(L::Liouvillian)::Liouvillian
  result = zero(L)
  for entry in term_pairs(L)
    action = first(entry)
    coefficient = last(entry)
    left_canonical = canonical_qadd(first(action))
    right_canonical = canonical_qadd(last(action))
    for left_entry in qadd_term_pairs(left_canonical)
      left_term = first(left_entry)
      left_coefficient = last(left_entry)
      for right_entry in qadd_term_pairs(right_canonical)
        right_term = first(right_entry)
        right_coefficient = last(right_entry)
        product = coefficient * left_coefficient * right_coefficient
        combined = simplify_coefficient(product)
        iszero(combined) && continue
        add_term!(
          result, monomial_operator(left_term), monomial_operator(right_term), combined
        )
      end
    end
  end
  return result
end

function multiply_coefficients(
  left::KossakowskiMatrix, right::KossakowskiMatrix
)::KossakowskiMatrix
  size(left, 2) == size(right, 1) ||
    throw(DimensionMismatch("matrix dimensions do not match"))
  result = coefficient_matrix(size(left, 1), size(right, 2))
  for row in axes(result, 1), column in axes(result, 2)
    value = coefficient_zero()
    for index in axes(left, 2)
      product = left[row, index] * right[index, column]
      value = value + product
    end
    result[row, column] = simplify_coefficient(value)
  end
  return result
end

function adjoint_coefficients(matrix::KossakowskiMatrix)::KossakowskiMatrix
  result = coefficient_matrix(size(matrix, 2), size(matrix, 1))
  for row in axes(matrix, 1), column in axes(matrix, 2)
    result[column, row] = conj(matrix[row, column])
  end
  return result
end

function sandwich_pivot_matrix_canonical(
  canonical::Liouvillian, frame::DissipativeFrame
)::KossakowskiMatrix
  q = length(frame.operators)
  result = coefficient_matrix(q, q)
  pivot_terms = frame.monomials[frame.pivot_rows]
  pivot_index = Dict{SQA.QTerm,Int}(
    term => index for (index, term) in enumerate(pivot_terms)
  )

  for entry in term_pairs(canonical)
    action = first(entry)
    coefficient = last(entry)
    left = first(action)
    right = last(action)
    (qadd_isone(left) || qadd_isone(right)) && continue
    left_term = first(first(qadd_term_pairs(left)))
    left_index = get(pivot_index, left_term, 0)
    iszero(left_index) && continue

    right_adjoint = canonical_qadd(adjoint(right))
    for right_entry in qadd_term_pairs(right_adjoint)
      right_term = first(right_entry)
      right_coefficient = last(right_entry)
      isempty(right_term.ops) && continue
      right_index = get(pivot_index, right_term, 0)
      iszero(right_index) && continue
      contribution = coefficient * conj(right_coefficient)
      result[left_index, right_index] = simplify_coefficient(
        result[left_index, right_index] + contribution
      )
    end
  end
  return result
end

function sandwich_pivot_matrix(L::Liouvillian, frame::DissipativeFrame)::KossakowskiMatrix
  return sandwich_pivot_matrix_canonical(canonical_liouvillian(L), frame)
end

function cross_dissipator(left::SQA.QAdd, right::SQA.QAdd)::Liouvillian
  identity = one(left)::SQA.QAdd
  right_adjoint = adjoint(right)::SQA.QAdd
  norm = (right_adjoint * left)::SQA.QAdd
  return (
    action(left, right_adjoint) +
    action(norm, identity, -1 // 2) +
    action(identity, norm, -1 // 2)
  )::Liouvillian
end

function dissipative_liouvillian(
  frame::DissipativeFrame, matrix::KossakowskiMatrix
)::Liouvillian
  q = length(frame.operators)
  size(matrix) == (q, q) ||
    throw(DimensionMismatch("Kossakowski matrix does not match frame"))
  result = zero(Liouvillian)
  for row in 1:q, column in 1:q
    coefficient = simplify_coefficient(matrix[row, column])
    iszero(coefficient) && continue
    dissipator_term = (
      coefficient * cross_dissipator(frame.operators[row], frame.operators[column])
    )::Liouvillian
    result = (result + dissipator_term)::Liouvillian
  end
  return canonical_liouvillian(result)
end

const FLOAT_ROUNDOFF_TOLERANCE = 1.0e-10

# Magnitude of a numeric constant and whether it is floating point; `nothing` if symbolic.
@inline function constant_magnitude(value)::Union{Nothing,Tuple{Float64,Bool}}
  value isa Int && return abs(Float64(value)), false
  value isa Rational{Int} && return abs(Float64(value)), false
  value isa Float64 && return abs(value), true
  value isa Complex{Int} && return abs(ComplexF64(value)), false
  value isa Complex{Rational{Int}} && return abs(ComplexF64(value)), false
  value isa ComplexF64 && return abs(value), true
  return nothing
end

function monomial_constant(term)::Tuple{Float64,Bool}
  constant = constant_magnitude(Symbolics.unwrap_const(term))
  constant === nothing || return constant
  Symbolics.iscall(term) || return 1.0, false
  # A monomial divided by a power of a symbol carries its prefactor in the numerator.
  Symbolics.operation(term) === (/) &&
    return monomial_constant(first(Symbolics.arguments(term)))
  Symbolics.operation(term) === (*) || return 1.0, false
  magnitude, inexact = 1.0, false
  for factor in Symbolics.arguments(term)
    factor_constant = constant_magnitude(Symbolics.unwrap_const(factor))
    factor_constant === nothing && continue
    magnitude *= first(factor_constant)
    inexact |= last(factor_constant)
  end
  return magnitude, inexact
end

function monomial_constants!(constants::Vector{Tuple{Float64,Bool}}, part::Symbolics.Num)
  expanded = Symbolics.unwrap(Symbolics.expand(part))
  terms = if Symbolics.iscall(expanded) && Symbolics.operation(expanded) === (+)
    Symbolics.arguments(expanded)
  else
    [expanded]
  end
  for term in terms
    constant = constant_magnitude(Symbolics.unwrap_const(term))
    constant !== nothing && iszero(first(constant)) && continue
    push!(constants, monomial_constant(term))
  end
  return constants
end

# Numeric prefactors of the monomials of a coefficient, with a flag for floating-point values.
function monomial_constants(value::SQA.CNum)::Vector{Tuple{Float64,Bool}}
  number = SQA.to_num(value)::Complex{Symbolics.Num}
  constants = Tuple{Float64,Bool}[]
  monomial_constants!(constants, real(number))
  monomial_constants!(constants, imag(number))
  return constants
end

# A difference is floating roundoff when every monomial carries a tiny floating prefactor.
function floating_roundoff(value::SQA.CNum, scale::Float64)::Bool
  constants = monomial_constants(value)
  isempty(constants) && return false
  return all(
    inexact && magnitude <= FLOAT_ROUNDOFF_TOLERANCE * scale for
    (magnitude, inexact) in constants
  )
end

function coefficient_scale(values)::Float64
  scale = 0.0
  for value in values
    for (magnitude, _) in monomial_constants(value)
      scale = max(scale, magnitude)
    end
  end
  return iszero(scale) ? 1.0 : scale
end

function hermitize_floating_roundoff!(matrix::KossakowskiMatrix)::KossakowskiMatrix
  scale = coefficient_scale(matrix)
  for row in axes(matrix, 1), column in row:size(matrix, 2)
    difference = simplify_coefficient(matrix[row, column] - conj(matrix[column, row]))
    (iszero(difference) || !floating_roundoff(difference, scale)) && continue
    average = simplify_coefficient((matrix[row, column] + conj(matrix[column, row])) / 2)
    matrix[row, column] = average
    matrix[column, row] = simplify_coefficient(conj(average))
  end
  return matrix
end

# Cheap sufficient test: the expanded difference vanishes term by term. This avoids the
# gcd-based `simplify`, which can overflow on large rational coefficients.
function expanded_conjugate_pair(left::SQA.CNum, right::SQA.CNum)::Bool
  difference =
    SQA.to_num(left)::Complex{Symbolics.Num} -
    conj(SQA.to_num(right)::Complex{Symbolics.Num})
  return iszero(Symbolics.expand(real(difference))) &&
         iszero(Symbolics.expand(imag(difference)))
end

function matrix_is_hermitian(matrix::KossakowskiMatrix)::Bool
  size(matrix, 1) == size(matrix, 2) || return false
  scale = coefficient_scale(matrix)
  for row in axes(matrix, 1), column in row:size(matrix, 2)
    expanded_conjugate_pair(matrix[row, column], matrix[column, row]) && continue
    difference = simplify_coefficient(matrix[row, column] - conj(matrix[column, row]))
    iszero(difference) || floating_roundoff(difference, scale) || return false
  end
  return true
end

function drop_floating_roundoff(residual::Liouvillian, scale::Float64)::Liouvillian
  result = zero(residual)
  for (left, right, coefficient) in terms(residual)
    floating_roundoff(coefficient, scale) && continue
    add_term!(result, left, right, coefficient)
  end
  return result
end

function canonical_has_two_sided_terms(canonical::Liouvillian)::Bool
  for entry in term_pairs(canonical)
    action = first(entry)
    left = first(action)
    right = last(action)
    (!qadd_isone(left) && !qadd_isone(right)) && return true
  end
  return false
end

function has_two_sided_terms(L::Liouvillian)::Bool
  return canonical_has_two_sided_terms(canonical_liouvillian(L))
end

function canonical_dissipative_support_terms(canonical::Liouvillian)::Vector{SQA.QTerm}
  terms_found = SQA.QTerm[]
  for entry in term_pairs(canonical)
    action = first(entry)
    left = first(action)
    right = last(action)
    (qadd_isone(left) || qadd_isone(right)) && continue
    left_term = first(first(qadd_term_pairs(left)))
    left_term in terms_found || push!(terms_found, left_term)
    right_adjoint = canonical_qadd(adjoint(right))
    for right_entry in qadd_term_pairs(right_adjoint)
      right_term = first(right_entry)
      isempty(right_term.ops) && continue
      right_term in terms_found || push!(terms_found, right_term)
    end
  end
  sort!(terms_found; by=SQA.term_order_key)
  return terms_found
end

function canonical_dissipative_support_directions(canonical::Liouvillian)::Vector{SQA.QAdd}
  terms_found = canonical_dissipative_support_terms(canonical)
  return [monomial_operator(term) for term in terms_found]
end

function dissipative_support_directions(L::Liouvillian)::Vector{SQA.QAdd}
  return canonical_dissipative_support_directions(canonical_liouvillian(L))
end

function residual_hamiltonian_canonical(canonical::Liouvillian)::SQA.QAdd
  H = zero(SQA.QAdd)
  for (left, right, coefficient) in terms(canonical)
    if !qadd_isone(left) && qadd_isone(right)
      H = (H + (im * coefficient) * left)::SQA.QAdd
    elseif qadd_isone(left) && !qadd_isone(right)
      continue
    else
      throw(GKSLCoordinateError("Liouvillian residual is not a Hamiltonian commutator"))
    end
  end
  H = SQA.simplify(H)::SQA.QAdd
  liouvillian_iszero(canonical_liouvillian(canonical - hamiltonian_action(H))) ||
    throw(GKSLCoordinateError("Liouvillian residual is not a Hamiltonian commutator"))
  hermiticity_residual = SQA.simplify(canonical_qadd(H - adjoint(H)))::SQA.QAdd
  qadd_iszero(hermiticity_residual) ||
    throw(GKSLCoordinateError("extracted Hamiltonian is not Hermitian modulo the identity"))
  return H
end

function residual_hamiltonian(residual::Liouvillian)::SQA.QAdd
  return residual_hamiltonian_canonical(canonical_liouvillian(residual))
end

function extract_gksl_canonical(
  canonical::Liouvillian, frame::DissipativeFrame
)::Tuple{SQA.QAdd,KossakowskiMatrix}
  sandwich = sandwich_pivot_matrix_canonical(canonical, frame)
  left = multiply_coefficients(frame.pivot_inverse, sandwich)
  matrix = multiply_coefficients(left, adjoint_coefficients(frame.pivot_inverse))
  simplify_matrix!(matrix)
  hermitize_floating_roundoff!(matrix)
  matrix_is_hermitian(matrix) || throw(
    GKSLCoordinateError(
      "extracted Kossakowski matrix is not Hermitian in the supplied frame"
    ),
  )

  dissipative = dissipative_liouvillian(frame, matrix)
  residual = canonical_liouvillian((canonical - dissipative)::Liouvillian)
  scale = coefficient_scale(coefficient for (_, _, coefficient) in terms(canonical))
  residual = drop_floating_roundoff(residual, scale)
  canonical_has_two_sided_terms(residual) && throw(
    ArgumentError("Liouvillian contains dissipative directions outside the supplied frame"),
  )
  H = residual_hamiltonian_canonical(residual)
  return H, matrix
end

function extract_gksl(
  L::Liouvillian, frame::DissipativeFrame
)::Tuple{SQA.QAdd,KossakowskiMatrix}
  return extract_gksl_canonical(canonical_liouvillian(L), frame)
end

function support_frame_canonical(canonical::Liouvillian)
  directions = canonical_dissipative_support_directions(canonical)
  return DissipativeFrame(Tuple(directions))
end

function support_frame(L::Liouvillian)
  return support_frame_canonical(canonical_liouvillian(L))
end

"""
    kossakowski(L::Liouvillian, frame::DissipativeFrame)

Return the Hermitian Kossakowski matrix of `L` in the ordered dissipative `frame`.

The matrix is extracted from the two-sided sandwich block and is exact in the symbolic SQA
algebra. The supplied frame must contain every dissipative direction of the Liouvillian.

See also [`DissipativeFrame`](@ref), [`kossakowski_component`](@ref), [`hamiltonian`](@ref).
"""
function kossakowski(L::Liouvillian, frame::DissipativeFrame)::KossakowskiMatrix
  _, matrix = extract_gksl(L, frame)
  return matrix
end

"""
    hamiltonian(L::Liouvillian[, frame::DissipativeFrame])

Return the coherent Hamiltonian of `L` in the canonical dissipative gauge, modulo an
additive multiple of the identity. With an explicit `frame`, the Liouvillian is
simultaneously checked to admit exact GKSL coordinates in that frame.

See also [`kossakowski`](@ref), [`hamiltonian_component`](@ref).
"""
function hamiltonian(L::Liouvillian, frame::DissipativeFrame)::SQA.QAdd
  H, _ = extract_gksl(L, frame)
  return H
end

function hamiltonian(L::Liouvillian)::SQA.QAdd
  canonical = canonical_liouvillian(L)
  canonical_has_two_sided_terms(canonical) ||
    return residual_hamiltonian_canonical(canonical)
  frame = support_frame_canonical(canonical)
  H, _ = extract_gksl_canonical(canonical, frame)
  return H
end

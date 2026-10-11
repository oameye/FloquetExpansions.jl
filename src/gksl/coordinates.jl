function canonical_liouvillian(L::Liouvillian)
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

function multiply_coefficients(left::KossakowskiMatrix, right::KossakowskiMatrix)
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

function adjoint_coefficients(matrix::KossakowskiMatrix)
  result = coefficient_matrix(size(matrix, 2), size(matrix, 1))
  for row in axes(matrix, 1), column in axes(matrix, 2)
    result[column, row] = conj(matrix[row, column])
  end
  return result
end

function sandwich_pivot_matrix_canonical(canonical::Liouvillian, frame::DissipativeFrame)
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

function sandwich_pivot_matrix(L::Liouvillian, frame::DissipativeFrame)
  return sandwich_pivot_matrix_canonical(canonical_liouvillian(L), frame)
end

function cross_dissipator(left::SQA.QAdd, right::SQA.QAdd)
  identity = one(left)
  right_adjoint = adjoint(right)
  norm = right_adjoint * left
  return action(left, right_adjoint) +
         action(norm, identity, -1 // 2) +
         action(identity, norm, -1 // 2)
end

function dissipative_liouvillian(frame::DissipativeFrame, matrix::KossakowskiMatrix)
  q = length(frame.operators)
  size(matrix) == (q, q) ||
    throw(DimensionMismatch("Kossakowski matrix does not match frame"))
  result = zero(Liouvillian)
  for row in 1:q, column in 1:q
    coefficient = simplify_coefficient(matrix[row, column])
    iszero(coefficient) && continue
    dissipator_term =
      coefficient * cross_dissipator(frame.operators[row], frame.operators[column])
    result = result + dissipator_term
  end
  return canonical_liouvillian(result)
end

const FLOAT_ROUNDOFF_TOLERANCE = 1.0e-10

const ExactConstant = Union{Int,Rational{Int},Complex{Int},Complex{Rational{Int}}}
const FloatConstant = Union{Float64,ComplexF64}
const NumericConstant = Union{ExactConstant,FloatConstant}

constant_magnitude(value::ExactConstant) = abs(ComplexF64(value)), false
constant_magnitude(value::FloatConstant) = abs(value), true

function monomial_constant(term)
  value = Symbolics.unwrap_const(term)
  value isa NumericConstant && return constant_magnitude(value)
  Symbolics.iscall(term) || return 1.0, false
  Symbolics.operation(term) === (/) &&
    return monomial_constant(first(Symbolics.arguments(term)))
  Symbolics.operation(term) === (*) || return 1.0, false
  magnitude, inexact = 1.0, false
  for factor in Symbolics.arguments(term)
    value = Symbolics.unwrap_const(factor)
    value isa NumericConstant || continue
    factor_magnitude, factor_inexact = constant_magnitude(value)
    magnitude *= factor_magnitude
    inexact |= factor_inexact
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
    value = Symbolics.unwrap_const(term)
    value isa NumericConstant && iszero(value) && continue
    push!(constants, monomial_constant(term))
  end
  return constants
end

function monomial_constants(value::SQA.CNum)
  number = SQA.to_num(value)
  constants = Tuple{Float64,Bool}[]
  monomial_constants!(constants, real(number))
  monomial_constants!(constants, imag(number))
  return constants
end

function has_floating_constant(node)
  Symbolics.unwrap_const(node) isa FloatConstant && return true
  Symbolics.iscall(node) || return false
  return any(has_floating_constant, Symbolics.arguments(node))
end

function has_floating_constant(value::SQA.CNum)
  number = SQA.to_num(value)
  return has_floating_constant(Symbolics.unwrap(real(number))) ||
         has_floating_constant(Symbolics.unwrap(imag(number)))
end

function floating_roundoff(value::SQA.CNum, reference, floor::Float64)
  has_floating_constant(value) || return false
  constants = monomial_constants(value)
  isempty(constants) && return false
  threshold = FLOAT_ROUNDOFF_TOLERANCE * max(coefficient_scale(reference), floor)
  for (magnitude, inexact) in constants
    (inexact && magnitude <= threshold) || return false
  end
  return true
end

function coefficient_scale(values)
  scale = 0.0
  for value in values
    for (magnitude, _) in monomial_constants(value)
      scale = max(scale, magnitude)
    end
  end
  return iszero(scale) ? 1.0 : scale
end

function hermitize_floating_roundoff!(matrix::KossakowskiMatrix)
  for row in axes(matrix, 1), column in row:size(matrix, 2)
    pair = (matrix[row, column], matrix[column, row])
    difference = simplify_coefficient(matrix[row, column] - conj(matrix[column, row]))
    (iszero(difference) || !floating_roundoff(difference, pair, 0.0)) && continue
    average = simplify_coefficient((matrix[row, column] + conj(matrix[column, row])) / 2)
    matrix[row, column] = average
    matrix[column, row] = simplify_coefficient(conj(average))
  end
  return matrix
end

function expanded_conjugate_pair(left::SQA.CNum, right::SQA.CNum)
  difference = SQA.to_num(left) - conj(SQA.to_num(right))
  return iszero(Symbolics.expand(real(difference))) &&
         iszero(Symbolics.expand(imag(difference)))
end

function matrix_is_hermitian(matrix::KossakowskiMatrix)
  size(matrix, 1) == size(matrix, 2) || return false
  for row in axes(matrix, 1), column in row:size(matrix, 2)
    expanded_conjugate_pair(matrix[row, column], matrix[column, row]) && continue
    pair = (matrix[row, column], matrix[column, row])
    difference = simplify_coefficient(matrix[row, column] - conj(matrix[column, row]))
    iszero(difference) || floating_roundoff(difference, pair, 0.0) || return false
  end
  return true
end

function drop_floating_roundoff(
  residual::Liouvillian,
  canonical::Liouvillian,
  dissipative::Liouvillian,
  matrix::KossakowskiMatrix,
)
  result = zero(residual)
  matrix_scale = any(has_floating_constant, matrix) ? coefficient_scale(matrix) : 0.0
  for (left, right, coefficient) in terms(residual)
    action = (left, right)
    reference = (
      get(canonical.terms, action, zero(SQA.CNum)),
      get(dissipative.terms, action, zero(SQA.CNum)),
    )
    two_sided = !qadd_isone(left) && !qadd_isone(right)
    floor = two_sided ? matrix_scale : 0.0
    floating_roundoff(coefficient, reference, floor) && continue
    add_term!(result, left, right, coefficient)
  end
  return result
end

function canonical_has_two_sided_terms(canonical::Liouvillian)
  for entry in term_pairs(canonical)
    action = first(entry)
    left = first(action)
    right = last(action)
    (!qadd_isone(left) && !qadd_isone(right)) && return true
  end
  return false
end

function has_two_sided_terms(L::Liouvillian)
  return canonical_has_two_sided_terms(canonical_liouvillian(L))
end

function canonical_dissipative_support_terms(canonical::Liouvillian)
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

function canonical_dissipative_support_directions(canonical::Liouvillian)
  terms_found = canonical_dissipative_support_terms(canonical)
  return [monomial_operator(term) for term in terms_found]
end

function dissipative_support_directions(L::Liouvillian)
  return canonical_dissipative_support_directions(canonical_liouvillian(L))
end

function residual_hamiltonian_canonical(canonical::Liouvillian)
  H = zero(SQA.QAdd)
  for (left, right, coefficient) in terms(canonical)
    if !qadd_isone(left) && qadd_isone(right)
      H = H + (im * coefficient) * left
    elseif qadd_isone(left) && !qadd_isone(right)
      continue
    else
      throw(GKSLCoordinateError("Liouvillian residual is not a Hamiltonian commutator"))
    end
  end
  H = SQA.simplify(H)
  liouvillian_iszero(canonical_liouvillian(canonical - hamiltonian_action(H))) ||
    throw(GKSLCoordinateError("Liouvillian residual is not a Hamiltonian commutator"))
  hermiticity_residual = SQA.simplify(canonical_qadd(H - adjoint(H)))
  qadd_iszero(hermiticity_residual) ||
    throw(GKSLCoordinateError("extracted Hamiltonian is not Hermitian modulo the identity"))
  return H
end

function residual_hamiltonian(residual::Liouvillian)
  return residual_hamiltonian_canonical(canonical_liouvillian(residual))
end

function extract_gksl_canonical(canonical::Liouvillian, frame::DissipativeFrame)
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
  residual = canonical_liouvillian(canonical - dissipative)
  residual = drop_floating_roundoff(residual, canonical, dissipative, matrix)
  canonical_has_two_sided_terms(residual) && throw(
    ArgumentError("Liouvillian contains dissipative directions outside the supplied frame"),
  )
  H = residual_hamiltonian_canonical(residual)
  return H, matrix
end

function extract_gksl(L::Liouvillian, frame::DissipativeFrame)
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
function kossakowski(L::Liouvillian, frame::DissipativeFrame)
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
function hamiltonian(L::Liouvillian, frame::DissipativeFrame)
  H, _ = extract_gksl(L, frame)
  return H
end

function hamiltonian(L::Liouvillian)
  canonical = canonical_liouvillian(L)
  canonical_has_two_sided_terms(canonical) ||
    return residual_hamiltonian_canonical(canonical)
  frame = support_frame_canonical(canonical)
  H, _ = extract_gksl_canonical(canonical, frame)
  return H
end

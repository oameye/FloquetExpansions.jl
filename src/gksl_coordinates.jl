const KossakowskiMatrix = Matrix{SQA.CNum}

"""
    DissipativeFrame(operators...)
    DissipativeFrame(operators::Tuple)

An ordered operator frame for Kossakowski coordinates. Operators are canonicalized with
SQA completeness relations and represented modulo the identity, since scalar shifts of a
jump direction change only the Hamiltonian gauge.

Ordering is part of the representation: two frames with the same span but different order
are distinct. The directions must be linearly independent modulo the identity and may be
nonorthogonal.

See also [`kossakowski`](@ref), [`hamiltonian`](@ref).
"""
struct DissipativeFrame
  operators::Vector{SQA.QAdd}
  monomials::Vector{SQA.QTerm}
  coordinates::KossakowskiMatrix
  pivot_rows::Vector{Int}
  pivot_inverse::KossakowskiMatrix
end

function Base.:(==)(left::DissipativeFrame, right::DissipativeFrame)
  return left.operators == right.operators
end
Base.isequal(left::DissipativeFrame, right::DissipativeFrame) =
  isequal(left.operators, right.operators)
Base.hash(frame::DissipativeFrame, h::UInt) =
  hash(:DissipativeFrame, hash(frame.operators, h))

@inline coefficient_zero() = convert(SQA.CNum, 0)
@inline coefficient_one() = convert(SQA.CNum, 1)
@inline simplify_coefficient(value::SQA.CNum) = SQA.simplify(value)

function canonical_qadd(operator::SQA.QField)
  result = SQA.simplify(SQA.expand_completeness(qadd(operator)))
  isempty(result.indices) || throw(
    ArgumentError("DissipativeFrame does not support operators with bound symbolic sums")
  )
  return result
end

function scalar_part(operator::SQA.QAdd)
  scalar = coefficient_zero()
  for (term, coefficient) in operator
    isempty(term.ops) || continue
    scalar = simplify_coefficient(scalar + coefficient)
  end
  return scalar
end

function projected_operator(operator::SQA.QField)
  canonical = canonical_qadd(operator)
  scalar = scalar_part(canonical)
  projected = iszero(scalar) ? canonical : SQA.simplify(canonical - scalar * one(canonical))
  iszero(projected) && throw(
    ArgumentError("a dissipative-frame direction cannot be proportional to the identity")
  )
  return projected
end

function sorted_term_pairs(operator::SQA.QAdd)
  pairs = collect(operator)
  sort!(pairs; by = pair -> SQA.term_order_key(first(pair)))
  return pairs
end

function monomial_operator(term::SQA.QTerm)
  arguments = SQA.QTermDict()
  arguments[term] = coefficient_one()
  return SQA.QAdd(arguments, SQA.Index[])
end

function coefficient_matrix(rows::Int, columns::Int)
  return fill(coefficient_zero(), rows, columns)
end

function simplify_matrix!(matrix::KossakowskiMatrix)
  for index in eachindex(matrix)
    matrix[index] = simplify_coefficient(matrix[index])
  end
  return matrix
end

function independent_pivot_rows(coordinates::KossakowskiMatrix)
  monomial_count, direction_count = size(coordinates)
  work = coefficient_matrix(direction_count, monomial_count)
  for row in 1:direction_count, column in 1:monomial_count
    work[row, column] = coordinates[column, row]
  end

  pivots = Int[]
  pivot_row = 1
  for column in 1:monomial_count
    candidate = 0
    for row in pivot_row:direction_count
      value = simplify_coefficient(work[row, column])
      work[row, column] = value
      if !iszero(value)
        candidate = row
        break
      end
    end
    iszero(candidate) && continue

    if candidate != pivot_row
      for trailing in 1:monomial_count
        work[pivot_row, trailing], work[candidate, trailing] =
          work[candidate, trailing], work[pivot_row, trailing]
      end
    end

    pivot = work[pivot_row, column]
    for row in (pivot_row + 1):direction_count
      entry = simplify_coefficient(work[row, column])
      iszero(entry) && continue
      factor = simplify_coefficient(entry * inv(pivot))
      for trailing in column:monomial_count
        work[row, trailing] = simplify_coefficient(
          work[row, trailing] - factor * work[pivot_row, trailing]
        )
      end
    end

    push!(pivots, column)
    pivot_row += 1
    pivot_row > direction_count && break
  end

  length(pivots) == direction_count || throw(
    ArgumentError("dissipative-frame directions are linearly dependent modulo the identity")
  )
  return pivots
end

function inverse_coefficients(matrix::KossakowskiMatrix)
  rows, columns = size(matrix)
  rows == columns || throw(ArgumentError("coefficient matrix must be square"))
  n = rows
  left = copy(matrix)
  right = coefficient_matrix(n, n)
  for index in 1:n
    right[index, index] = coefficient_one()
  end

  for column in 1:n
    candidate = 0
    for row in column:n
      value = simplify_coefficient(left[row, column])
      left[row, column] = value
      if !iszero(value)
        candidate = row
        break
      end
    end
    iszero(candidate) && throw(
      ArgumentError("dissipative-frame coordinate pivot is structurally singular")
    )

    if candidate != column
      for trailing in 1:n
        left[column, trailing], left[candidate, trailing] =
          left[candidate, trailing], left[column, trailing]
        right[column, trailing], right[candidate, trailing] =
          right[candidate, trailing], right[column, trailing]
      end
    end

    pivot_inverse = inv(left[column, column])
    for trailing in 1:n
      left[column, trailing] = simplify_coefficient(left[column, trailing] * pivot_inverse)
      right[column, trailing] = simplify_coefficient(right[column, trailing] * pivot_inverse)
    end

    for row in 1:n
      row == column && continue
      factor = simplify_coefficient(left[row, column])
      iszero(factor) && continue
      for trailing in 1:n
        left[row, trailing] = simplify_coefficient(
          left[row, trailing] - factor * left[column, trailing]
        )
        right[row, trailing] = simplify_coefficient(
          right[row, trailing] - factor * right[column, trailing]
        )
      end
    end
  end
  return right
end

function build_dissipative_frame(operators)
  isempty(operators) && throw(ArgumentError("DissipativeFrame requires at least one direction"))
  projected = SQA.QAdd[]
  sizehint!(projected, length(operators))
  for operator in operators
    operator isa SQA.QField || throw(
      ArgumentError("every dissipative-frame direction must be an SQA operator expression")
    )
    push!(projected, projected_operator(operator))
  end

  monomials = SQA.QTerm[]
  for operator in projected
    for (term, _) in sorted_term_pairs(operator)
      isempty(term.ops) && continue
      term in monomials || push!(monomials, term)
    end
  end

  coordinates = coefficient_matrix(length(monomials), length(projected))
  row_index = Dict{SQA.QTerm,Int}(term => row for (row, term) in enumerate(monomials))
  for (column, operator) in enumerate(projected), (term, coefficient) in operator
    isempty(term.ops) && continue
    coordinates[row_index[term], column] = coefficient
  end

  pivots = independent_pivot_rows(coordinates)
  pivot_matrix = coefficient_matrix(length(pivots), length(projected))
  for row in eachindex(pivots), column in eachindex(projected)
    pivot_matrix[row, column] = coordinates[pivots[row], column]
  end
  pivot_inverse = inverse_coefficients(pivot_matrix)

  return DissipativeFrame(projected, monomials, coordinates, pivots, pivot_inverse)
end

DissipativeFrame(operators::Tuple) = build_dissipative_frame(operators)
DissipativeFrame(operators::SQA.QField...) = build_dissipative_frame(operators)

function canonical_liouvillian(L::Liouvillian)
  result = zero(L)
  for (left, right, coefficient) in terms(L)
    left_canonical = canonical_qadd(left)
    right_canonical = canonical_qadd(right)
    for (left_term, left_coefficient) in left_canonical,
      (right_term, right_coefficient) in right_canonical

      combined = simplify_coefficient(coefficient * left_coefficient * right_coefficient)
      iszero(combined) && continue
      add_term!(
        result,
        monomial_operator(left_term),
        monomial_operator(right_term),
        combined,
      )
    end
  end
  return result
end

function multiply_coefficients(
  left::KossakowskiMatrix, right::KossakowskiMatrix
)
  size(left, 2) == size(right, 1) || throw(DimensionMismatch("matrix dimensions do not match"))
  result = coefficient_matrix(size(left, 1), size(right, 2))
  for row in axes(result, 1), column in axes(result, 2)
    value = coefficient_zero()
    for index in axes(left, 2)
      value = value + left[row, index] * right[index, column]
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

function sandwich_pivot_matrix(L::Liouvillian, frame::DissipativeFrame)
  q = length(frame.operators)
  result = coefficient_matrix(q, q)
  pivot_terms = frame.monomials[frame.pivot_rows]
  pivot_index = Dict{SQA.QTerm,Int}(term => index for (index, term) in enumerate(pivot_terms))

  for (left, right, coefficient) in terms(canonical_liouvillian(L))
    (isone(left) || isone(right)) && continue
    left_term = first(first(left))
    left_index = get(pivot_index, left_term, 0)
    iszero(left_index) && continue

    right_adjoint = canonical_qadd(adjoint(right))
    for (right_term, right_coefficient) in right_adjoint
      isempty(right_term.ops) && continue
      right_index = get(pivot_index, right_term, 0)
      iszero(right_index) && continue
      result[left_index, right_index] = simplify_coefficient(
        result[left_index, right_index] + coefficient * conj(right_coefficient)
      )
    end
  end
  return result
end

function cross_dissipator(left::SQA.QAdd, right::SQA.QAdd)
  identity = one(left)
  norm = adjoint(right) * left
  return action(left, adjoint(right)) +
         action(norm, identity, -1 // 2) +
         action(identity, norm, -1 // 2)
end

function dissipative_liouvillian(frame::DissipativeFrame, matrix::KossakowskiMatrix)
  q = length(frame.operators)
  size(matrix) == (q, q) || throw(DimensionMismatch("Kossakowski matrix does not match frame"))
  result = zero(Liouvillian)
  for row in 1:q, column in 1:q
    coefficient = simplify_coefficient(matrix[row, column])
    iszero(coefficient) && continue
    result = result + coefficient * cross_dissipator(frame.operators[row], frame.operators[column])
  end
  return canonical_liouvillian(result)
end

function matrix_is_hermitian(matrix::KossakowskiMatrix)
  size(matrix, 1) == size(matrix, 2) || return false
  for row in axes(matrix, 1), column in row:size(matrix, 2)
    difference = simplify_coefficient(matrix[row, column] - conj(matrix[column, row]))
    iszero(difference) || return false
  end
  return true
end

function has_two_sided_terms(L::Liouvillian)
  for (left, right, _) in terms(canonical_liouvillian(L))
    (!isone(left) && !isone(right)) && return true
  end
  return false
end

function residual_hamiltonian(residual::Liouvillian)
  canonical = canonical_liouvillian(residual)
  H = zero(SQA.QAdd)
  for (left, right, coefficient) in terms(canonical)
    if !isone(left) && isone(right)
      H = H + (im * coefficient) * left
    elseif isone(left) && !isone(right)
      continue
    else
      throw(ArgumentError("Liouvillian residual is not a Hamiltonian commutator"))
    end
  end
  H = SQA.simplify(H)
  iszero(canonical_liouvillian(canonical - hamiltonian_action(H))) || throw(
    ArgumentError("Liouvillian residual is not a Hamiltonian commutator")
  )
  iszero(SQA.simplify(canonical_qadd(H - adjoint(H)))) || throw(
    ArgumentError("extracted Hamiltonian is not Hermitian modulo the identity")
  )
  return H
end

function extract_gksl(L::Liouvillian, frame::DissipativeFrame)
  sandwich = sandwich_pivot_matrix(L, frame)
  left = multiply_coefficients(frame.pivot_inverse, sandwich)
  matrix = multiply_coefficients(left, adjoint_coefficients(frame.pivot_inverse))
  simplify_matrix!(matrix)
  matrix_is_hermitian(matrix) || throw(
    ArgumentError("extracted Kossakowski matrix is not Hermitian in the supplied frame")
  )

  dissipative = dissipative_liouvillian(frame, matrix)
  residual = canonical_liouvillian(L - dissipative)
  has_two_sided_terms(residual) && throw(
    ArgumentError("Liouvillian contains dissipative directions outside the supplied frame")
  )
  H = residual_hamiltonian(residual)
  return H, matrix
end

function support_frame(L::Liouvillian)
  canonical = canonical_liouvillian(L)
  terms_found = SQA.QTerm[]
  for (left, right, _) in terms(canonical)
    (isone(left) || isone(right)) && continue
    left_term = first(first(left))
    left_term in terms_found || push!(terms_found, left_term)
    for (right_term, _) in canonical_qadd(adjoint(right))
      isempty(right_term.ops) && continue
      right_term in terms_found || push!(terms_found, right_term)
    end
  end
  sort!(terms_found; by = SQA.term_order_key)
  return DissipativeFrame(Tuple(monomial_operator(term) for term in terms_found))
end

"""
    kossakowski(L::Liouvillian, frame::DissipativeFrame)
    kossakowski(expansion::FloquetExpansion, frame::DissipativeFrame)

Return the Hermitian Kossakowski matrix of a Liouvillian, or of the finite effective
Liouvillian of a Floquet expansion, in the ordered dissipative `frame`.

The matrix is extracted from the two-sided sandwich block and is exact in the symbolic SQA
algebra. The supplied frame must contain every dissipative direction of the Liouvillian.

See also [`DissipativeFrame`](@ref), [`kossakowski_component`](@ref), [`hamiltonian`](@ref).
"""
function kossakowski(L::Liouvillian, frame::DissipativeFrame)
  _, matrix = extract_gksl(L, frame)
  return matrix
end

function kossakowski(
  expansion::FloquetExpansion{G,P,E}, frame::DissipativeFrame
) where {G,P,E<:Liouvillian}
  return kossakowski(effective_generator(expansion), frame)
end

"""
    kossakowski_component(expansion::FloquetExpansion, frame::DissipativeFrame, n::Int)

Return the order-`n` Kossakowski contribution of a Liouvillian Floquet expansion in the
ordered dissipative `frame`, including the corresponding inverse-drive-frequency scaling.

See also [`kossakowski`](@ref), [`effective_generator`](@ref).
"""
function kossakowski_component(
  expansion::FloquetExpansion{G,P,E}, frame::DissipativeFrame, n::Int
) where {G,P,E<:Liouvillian}
  return kossakowski(effective_generator(expansion, n), frame)
end

"""
    hamiltonian(L::Liouvillian[, frame::DissipativeFrame])
    hamiltonian(expansion::FloquetExpansion)

Return the coherent Hamiltonian in the canonical dissipative gauge, modulo an additive
multiple of the identity. With an explicit `frame`, the Liouvillian is simultaneously
checked to admit exact GKSL coordinates in that frame.

For Hamiltonian Floquet expansions this is the effective Hamiltonian itself. For Liouvillian
expansions the dissipative sandwich support is removed before extracting the commutator.

See also [`hamiltonian_component`](@ref), [`kossakowski`](@ref).
"""
function hamiltonian(L::Liouvillian, frame::DissipativeFrame)
  H, _ = extract_gksl(L, frame)
  return H
end

function hamiltonian(L::Liouvillian)
  canonical = canonical_liouvillian(L)
  has_two_sided_terms(canonical) || return residual_hamiltonian(canonical)
  frame = support_frame(canonical)
  return hamiltonian(canonical, frame)
end

hamiltonian(expansion::FloquetExpansion{G,P,E}) where {G,P,E<:SQA.QAdd} =
  effective_generator(expansion)
hamiltonian(expansion::FloquetExpansion{G,P,E}) where {G,P,E<:Liouvillian} =
  hamiltonian(effective_generator(expansion))

"""
    hamiltonian_component(expansion::FloquetExpansion, n::Int)

Return the coherent order-`n` contribution of a Floquet expansion, including its
inverse-drive-frequency scaling. The Hamiltonian is defined modulo an additive scalar
multiple of the identity.

See also [`hamiltonian`](@ref), [`effective_generator`](@ref).
"""
hamiltonian_component(
  expansion::FloquetExpansion{G,P,E}, n::Int
) where {G,P,E<:SQA.QAdd} = effective_generator(expansion, n)

function hamiltonian_component(
  expansion::FloquetExpansion{G,P,E}, n::Int
) where {G,P,E<:Liouvillian}
  return hamiltonian(effective_generator(expansion, n))
end

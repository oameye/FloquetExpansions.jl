from pathlib import Path


def replace_once(text, old, new, label):
    if text.count(old) != 1:
        raise RuntimeError(f"{label}: expected exactly one match, found {text.count(old)}")
    return text.replace(old, new, 1)


# ---------------------------------------------------------------------------
# GKSL coordinates: share coordinate/support construction and keep canonical
# inputs canonical through extraction instead of simplifying them repeatedly.
# ---------------------------------------------------------------------------
path = Path("src/gksl_coordinates.jl")
text = path.read_text()

old = '''function build_dissipative_frame(operators)
  isempty(operators) &&
    throw(ArgumentError("DissipativeFrame requires at least one direction"))
  projected = frame_operator_tuple(operators)

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
'''
new = '''function frame_coordinate_data(operators)
  monomials = SQA.QTerm[]
  for operator in operators
    for (term, _) in sorted_term_pairs(operator)
      isempty(term.ops) && continue
      term in monomials || push!(monomials, term)
    end
  end

  coordinates = coefficient_matrix(length(monomials), length(operators))
  row_index = Dict{SQA.QTerm,Int}(term => row for (row, term) in enumerate(monomials))
  for (column, operator) in enumerate(operators), (term, coefficient) in operator
    isempty(term.ops) && continue
    coordinates[row_index[term], column] = coefficient
  end
  return monomials, coordinates
end

function frame_direction_coordinates(operators::Vector{SQA.QAdd})
  _, coordinates = frame_coordinate_data(operators)
  return coordinates
end

function build_dissipative_frame(operators)
  isempty(operators) &&
    throw(ArgumentError("DissipativeFrame requires at least one direction"))
  projected = frame_operator_tuple(operators)
  monomials, coordinates = frame_coordinate_data(projected)

  pivots = independent_pivot_rows(coordinates)
'''
text = replace_once(text, old, new, "frame coordinate data")

old = '''function sandwich_pivot_matrix(L::Liouvillian, frame::DissipativeFrame)::KossakowskiMatrix
  q = length(frame.operators)
  result = coefficient_matrix(q, q)
  pivot_terms = frame.monomials[frame.pivot_rows]
  pivot_index = Dict{SQA.QTerm,Int}(
    term => index for (index, term) in enumerate(pivot_terms)
  )

  canonical = canonical_liouvillian(L)
  for entry in term_pairs(canonical)
'''
new = '''function sandwich_pivot_matrix_canonical(
  canonical::Liouvillian, frame::DissipativeFrame
)::KossakowskiMatrix
  q = length(frame.operators)
  result = coefficient_matrix(q, q)
  pivot_terms = frame.monomials[frame.pivot_rows]
  pivot_index = Dict{SQA.QTerm,Int}(
    term => index for (index, term) in enumerate(pivot_terms)
  )

  for entry in term_pairs(canonical)
'''
text = replace_once(text, old, new, "canonical sandwich header")

old = '''  return result
end

function cross_dissipator(left::SQA.QAdd, right::SQA.QAdd)::Liouvillian
'''
new = '''  return result
end

function sandwich_pivot_matrix(L::Liouvillian, frame::DissipativeFrame)::KossakowskiMatrix
  return sandwich_pivot_matrix_canonical(canonical_liouvillian(L), frame)
end

function cross_dissipator(left::SQA.QAdd, right::SQA.QAdd)::Liouvillian
'''
text = replace_once(text, old, new, "canonical sandwich wrapper")

old = '''function has_two_sided_terms(L::Liouvillian)::Bool
  for (left, right, _) in terms(canonical_liouvillian(L))
    (!qadd_isone(left) && !qadd_isone(right)) && return true
  end
  return false
end

function residual_hamiltonian(residual::Liouvillian)::SQA.QAdd
  canonical = canonical_liouvillian(residual)
  H = zero(SQA.QAdd)
'''
new = '''function canonical_has_two_sided_terms(canonical::Liouvillian)::Bool
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
'''
text = replace_once(text, old, new, "canonical support and residual")

old = '''  return H
end

function extract_gksl(
  L::Liouvillian, frame::DissipativeFrame
)::Tuple{SQA.QAdd,KossakowskiMatrix}
  sandwich = sandwich_pivot_matrix(L, frame)
'''
new = '''  return H
end

function residual_hamiltonian(residual::Liouvillian)::SQA.QAdd
  return residual_hamiltonian_canonical(canonical_liouvillian(residual))
end

function extract_gksl_canonical(
  canonical::Liouvillian, frame::DissipativeFrame
)::Tuple{SQA.QAdd,KossakowskiMatrix}
  sandwich = sandwich_pivot_matrix_canonical(canonical, frame)
'''
text = replace_once(text, old, new, "canonical extract header")

old = '''  dissipative = dissipative_liouvillian(frame, matrix)
  residual = canonical_liouvillian((L - dissipative)::Liouvillian)
  has_two_sided_terms(residual) && throw(
    ArgumentError("Liouvillian contains dissipative directions outside the supplied frame"),
  )
  H = residual_hamiltonian(residual)
  return H, matrix
end

function support_frame(L::Liouvillian)
  canonical = canonical_liouvillian(L)
  terms_found = SQA.QTerm[]
  for (left, right, _) in terms(canonical)
    (qadd_isone(left) || qadd_isone(right)) && continue
    left_term = first(first(left))
    left_term in terms_found || push!(terms_found, left_term)
    for (right_term, _) in canonical_qadd(adjoint(right))
      isempty(right_term.ops) && continue
      right_term in terms_found || push!(terms_found, right_term)
    end
  end
  sort!(terms_found; by=SQA.term_order_key)
  return DissipativeFrame(Tuple(monomial_operator(term) for term in terms_found))
end
'''
new = '''  dissipative = dissipative_liouvillian(frame, matrix)
  residual = canonical_liouvillian((canonical - dissipative)::Liouvillian)
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
'''
text = replace_once(text, old, new, "canonical extract residual and support")

old = '''function hamiltonian(L::Liouvillian)::SQA.QAdd
  canonical = canonical_liouvillian(L)
  has_two_sided_terms(canonical) || return residual_hamiltonian(canonical)
  frame = support_frame(canonical)
  return hamiltonian(canonical, frame)
end
'''
new = '''function hamiltonian(L::Liouvillian)::SQA.QAdd
  canonical = canonical_liouvillian(L)
  canonical_has_two_sided_terms(canonical) || return residual_hamiltonian_canonical(canonical)
  frame = support_frame_canonical(canonical)
  H, _ = extract_gksl_canonical(canonical, frame)
  return H
end
'''
text = replace_once(text, old, new, "hamiltonian canonical path")
path.write_text(text)


# ---------------------------------------------------------------------------
# Automatic frame discovery: use the shared coordinate/support kernels and do
# not canonicalize once for a pre-scan and again for support extraction.
# ---------------------------------------------------------------------------
path = Path("src/completion_frame.jl")
text = path.read_text()
old = '''function frame_direction_coordinates(operators::Vector{SQA.QAdd})
  monomials = SQA.QTerm[]
  for operator in operators
    for (term, _) in sorted_term_pairs(operator)
      isempty(term.ops) && continue
      term in monomials || push!(monomials, term)
    end
  end

  coordinates = coefficient_matrix(length(monomials), length(operators))
  row_index = Dict{SQA.QTerm,Int}(term => row for (row, term) in enumerate(monomials))
  for (column, operator) in enumerate(operators), (term, coefficient) in operator
    isempty(term.ops) && continue
    coordinates[row_index[term], column] = coefficient
  end
  return coordinates
end

'''
text = replace_once(text, old, "", "remove duplicate frame coordinates")
old = '''function generated_support_directions(L::Liouvillian)
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
  sort!(terms_found; by=SQA.term_order_key)
  return [monomial_operator(term) for term in terms_found]
end

function append_generated_directions!(operators::Vector{SQA.QAdd}, L::Liouvillian)
  has_two_sided_terms(L) || return operators
  for operator in generated_support_directions(L)
    append_frame_candidate!(operators, operator)
  end
  return operators
end
'''
new = '''function append_generated_directions!(operators::Vector{SQA.QAdd}, L::Liouvillian)
  for operator in dissipative_support_directions(L)
    append_frame_candidate!(operators, operator)
  end
  return operators
end
'''
text = replace_once(text, old, new, "share dissipative support")
path.write_text(text)


# ---------------------------------------------------------------------------
# Matrix-series core: avoid captured generator closures, redundant matrix
# allocation, and encode the LDL matrix coordinate as one semantic argument.
# ---------------------------------------------------------------------------
path = Path("src/matrix_series.jl")
text = path.read_text()
old = '''function condition_contains(conditions::Vector{CompletionScalar}, x::CompletionScalar)
  return any(p -> structurally_equal(p, x), conditions)
end
'''
new = '''function condition_contains(conditions::Vector{CompletionScalar}, x::CompletionScalar)
  for condition in conditions
    structurally_equal(condition, x) && return true
  end
  return false
end
'''
text = replace_once(text, old, new, "condition lookup loop")
old = '''function validate_matrix_series(series::MatrixSeries)
  isempty(series) &&
    throw(ArgumentError("matrix series must contain at least one coefficient"))
  dims = size(first(series))
  all(size(A) == dims for A in series) ||
    throw(DimensionMismatch("all matrix-series coefficients must have the same dimensions"))
  return dims
end
'''
new = '''function validate_matrix_series(series::MatrixSeries)
  isempty(series) &&
    throw(ArgumentError("matrix series must contain at least one coefficient"))
  dims = size(first(series))
  for matrix in series
    size(matrix) == dims ||
      throw(DimensionMismatch("all matrix-series coefficients must have the same dimensions"))
  end
  return dims
end
'''
text = replace_once(text, old, new, "matrix series validation loop")
old = '''  result = [completion_matrix_zeros(m, p) for _ in 0:N]
  for n in 0:N
    coefficient = completion_matrix_zeros(m, p)
'''
new = '''  result = Vector{CompletionMatrix}(undef, N + 1)
  for n in 0:N
    coefficient = completion_matrix_zeros(m, p)
'''
text = replace_once(text, old, new, "matrix series result allocation")
old = '''function graded_ldl_numerator(
  A::MatrixSeries,
  lower::MatrixSeries,
  diagonal::Vector{ScalarSeries},
  i::Int,
  j::Int,
  N::Int,
)::ScalarSeries
  numerator = scalar_series_entry(A, i, j, N)
'''
new = '''function graded_ldl_numerator(
  A::MatrixSeries,
  lower::MatrixSeries,
  diagonal::Vector{ScalarSeries},
  position::CartesianIndex{2},
  N::Int,
)::ScalarSeries
  i, j = Tuple(position)
  numerator = scalar_series_entry(A, i, j, N)
'''
text = replace_once(text, old, new, "LDL coordinate argument")
old = '''      numerator = graded_ldl_numerator(A, lower, diagonal, i, j, N)
'''
new = '''      numerator = graded_ldl_numerator(A, lower, diagonal, CartesianIndex(i, j), N)
'''
text = replace_once(text, old, new, "LDL coordinate call")
path.write_text(text)

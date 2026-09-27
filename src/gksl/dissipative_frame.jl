const KossakowskiMatrix = Matrix{SQA.CNum}

"""
    GKSLCoordinateError

Raised when exact GKSL coordinates cannot be constructed in a requested dissipative frame.
"""
struct GKSLCoordinateError <: Exception
  message::String
end

Base.showerror(io::IO, error::GKSLCoordinateError) = print(io, error.message)

"""
    DissipativeFrame(operators...)
    DissipativeFrame(operators::Tuple)
    DissipativeFrame(operators::AbstractVector)

An ordered operator frame for Kossakowski coordinates. Operators are canonicalized with
SQA completeness relations and represented modulo the identity, since scalar shifts of a
jump direction change only the Hamiltonian gauge.

Ordering is part of the representation: two frames with the same span but different order
are distinct. The directions must be linearly independent modulo the identity and may be
nonorthogonal.

See also [`kossakowski`](@ref), [`hamiltonian`](@ref).
"""
struct DissipativeFrame{O<:Tuple}
  operators::O
  monomials::Vector{SQA.QTerm}
  coordinates::KossakowskiMatrix
  pivot_rows::Vector{Int}
  pivot_inverse::KossakowskiMatrix
end

function Base.:(==)(left::DissipativeFrame, right::DissipativeFrame)
  return left.operators == right.operators
end
function Base.isequal(left::DissipativeFrame, right::DissipativeFrame)
  return isequal(left.operators, right.operators)
end
function Base.hash(frame::DissipativeFrame, h::UInt)
  return hash(:DissipativeFrame, hash(frame.operators, h))
end

function Base.show(io::IO, frame::DissipativeFrame)
  print(io, "DissipativeFrame(")
  for (index, operator) in enumerate(frame.operators)
    index == 1 || print(io, ", ")
    show(io, operator)
  end
  return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", frame::DissipativeFrame)
  count = length(frame.operators)
  label = count == 1 ? "direction" : "directions"
  print(io, "DissipativeFrame with ", count, " ordered ", label, ":")
  for (index, operator) in enumerate(frame.operators)
    print(io, "\n  ", index, ": ")
    show(io, operator)
  end
  return nothing
end

@inline coefficient_zero()::SQA.CNum = convert(SQA.CNum, 0)
@inline coefficient_one()::SQA.CNum = convert(SQA.CNum, 1)
@inline function simplify_coefficient(value::SQA.CNum)::SQA.CNum
  return SQA.simplify(value)::SQA.CNum
end
@inline qadd_isone(value::SQA.QAdd)::Bool = isone(value)::Bool
@inline qadd_iszero(value::SQA.QAdd)::Bool = iszero(value)::Bool
@inline liouvillian_iszero(value::Liouvillian)::Bool = iszero(value)::Bool
@inline qadd_term_pairs(value::SQA.QAdd) = pairs(getfield(value, :arguments))

function canonical_qadd(operator::SQA.QField)::SQA.QAdd
  result = SQA.simplify(SQA.expand_completeness(qadd(operator)))::SQA.QAdd
  isempty(result.indices) || throw(
    ArgumentError("GKSL coordinates do not support operators with bound symbolic sums")
  )
  return result
end

function scalar_part(operator::SQA.QAdd)::SQA.CNum
  scalar = coefficient_zero()
  for (term, coefficient) in operator
    isempty(term.ops) || continue
    scalar = simplify_coefficient(scalar + coefficient)
  end
  return scalar
end

function projected_operator(operator::SQA.QField)::SQA.QAdd
  canonical = canonical_qadd(operator)
  scalar = scalar_part(canonical)
  projected = if iszero(scalar)
    canonical
  else
    SQA.simplify(canonical - scalar * one(canonical))::SQA.QAdd
  end
  qadd_iszero(projected) && throw(
    ArgumentError("a dissipative-frame direction cannot be proportional to the identity")
  )
  return projected
end

project_frame_operator(operator::SQA.QField)::SQA.QAdd = projected_operator(operator)
function project_frame_operator(::Any)
  return throw(
    ArgumentError("every dissipative-frame direction must be an SQA operator expression")
  )
end

function frame_operator_tuple(operators::Tuple{Vararg{SQA.QField,N}}) where {N}
  return ntuple(index -> project_frame_operator(operators[index]), Val(N))
end
frame_operator_tuple(operators::Tuple) = map(project_frame_operator, operators)
function frame_operator_tuple(operators::AbstractVector)
  return Tuple(project_frame_operator(operator) for operator in operators)
end

function sorted_term_pairs(operator::SQA.QAdd)
  pairs = collect(operator)
  sort!(pairs; by=pair -> SQA.term_order_key(first(pair)))
  return pairs
end

function monomial_operator(term::SQA.QTerm)::SQA.QAdd
  arguments = SQA.QTermDict()
  arguments[term] = coefficient_one()
  return SQA.QAdd(arguments, SQA.Index[])
end

function coefficient_matrix(rows::Int, columns::Int)::KossakowskiMatrix
  return fill(coefficient_zero(), rows, columns)
end

function simplify_matrix!(matrix::KossakowskiMatrix)::KossakowskiMatrix
  for index in eachindex(matrix)
    matrix[index] = simplify_coefficient(matrix[index])
  end
  return matrix
end

function transposed_coordinate_matrix(coordinates::KossakowskiMatrix)
  monomial_count, direction_count = size(coordinates)
  work = coefficient_matrix(direction_count, monomial_count)
  for row in 1:direction_count, column in 1:monomial_count
    work[row, column] = coordinates[column, row]
  end
  return work
end

function coefficient_pivot_row!(work::KossakowskiMatrix, pivot::CartesianIndex{2})::Int
  first_row, column = Tuple(pivot)
  for row in first_row:size(work, 1)
    value = simplify_coefficient(work[row, column])
    work[row, column] = value
    iszero(value) || return row
  end
  return 0
end

function swap_coefficient_rows!(matrix::KossakowskiMatrix, rows::NTuple{2,Int})
  first_row, second_row = rows
  first_row == second_row && return matrix
  for column in axes(matrix, 2)
    matrix[first_row, column], matrix[second_row, column] = matrix[second_row, column],
    matrix[first_row, column]
  end
  return matrix
end

function eliminate_coordinate_column!(
  work::KossakowskiMatrix, pivot_index::CartesianIndex{2}
)
  pivot_row, column = Tuple(pivot_index)
  direction_count, monomial_count = size(work)
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
  return work
end

function coordinate_pivot_rows(coordinates::KossakowskiMatrix)::Vector{Int}
  monomial_count, direction_count = size(coordinates)
  direction_count == 0 && return Int[]
  monomial_count < direction_count && return Int[]

  work = transposed_coordinate_matrix(coordinates)
  pivots = Int[]
  pivot_row = 1
  for column in 1:monomial_count
    pivot_index = CartesianIndex(pivot_row, column)
    candidate = coefficient_pivot_row!(work, pivot_index)
    iszero(candidate) && continue
    swap_coefficient_rows!(work, (pivot_row, candidate))
    eliminate_coordinate_column!(work, pivot_index)
    push!(pivots, column)
    pivot_row += 1
    pivot_row > direction_count && break
  end
  return pivots
end

function independent_pivot_rows(coordinates::KossakowskiMatrix)
  _, direction_count = size(coordinates)
  pivots = coordinate_pivot_rows(coordinates)
  length(pivots) == direction_count || throw(
    ArgumentError(
      "dissipative-frame directions are linearly dependent modulo the identity"
    ),
  )
  return pivots
end

function coefficient_identity(n::Int)::KossakowskiMatrix
  result = coefficient_matrix(n, n)
  for index in 1:n
    result[index, index] = coefficient_one()
  end
  return result
end

struct CoefficientInverseWorkspace
  left::KossakowskiMatrix
  right::KossakowskiMatrix
end

function normalize_inverse_pivot_row!(workspace::CoefficientInverseWorkspace, column::Int)
  left = workspace.left
  right = workspace.right
  pivot_inverse = inv(left[column, column])
  for trailing in axes(left, 2)
    left[column, trailing] = simplify_coefficient(left[column, trailing] * pivot_inverse)
    right[column, trailing] = simplify_coefficient(right[column, trailing] * pivot_inverse)
  end
  return nothing
end

function eliminate_inverse_column!(workspace::CoefficientInverseWorkspace, column::Int)
  left = workspace.left
  right = workspace.right
  for row in axes(left, 1)
    row == column && continue
    factor = simplify_coefficient(left[row, column])
    iszero(factor) && continue
    for trailing in axes(left, 2)
      left[row, trailing] = simplify_coefficient(
        left[row, trailing] - factor * left[column, trailing]
      )
      right[row, trailing] = simplify_coefficient(
        right[row, trailing] - factor * right[column, trailing]
      )
    end
  end
  return nothing
end

function inverse_coefficients(matrix::KossakowskiMatrix)::KossakowskiMatrix
  rows, columns = size(matrix)
  rows == columns || throw(DimensionMismatch("coefficient matrix must be square"))
  n = rows
  left = copy(matrix)
  right = coefficient_identity(n)
  workspace = CoefficientInverseWorkspace(left, right)

  for column in 1:n
    candidate = coefficient_pivot_row!(left, CartesianIndex(column, column))
    iszero(candidate) && throw(
      GKSLCoordinateError("dissipative-frame coordinate pivot is structurally singular")
    )
    swap_coefficient_rows!(left, (column, candidate))
    swap_coefficient_rows!(right, (column, candidate))
    normalize_inverse_pivot_row!(workspace, column)
    eliminate_inverse_column!(workspace, column)
  end
  return right
end

function frame_coordinate_data(operators)
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
  return assemble_dissipative_frame(operators)
end

function assemble_dissipative_frame(operators)
  projected = frame_operator_tuple(operators)
  monomials, coordinates = frame_coordinate_data(projected)

  pivots = independent_pivot_rows(coordinates)
  pivot_matrix = coefficient_matrix(length(pivots), length(projected))
  for row in eachindex(pivots), column in eachindex(projected)
    pivot_matrix[row, column] = coordinates[pivots[row], column]
  end
  pivot_inverse = inverse_coefficients(pivot_matrix)

  return DissipativeFrame(projected, monomials, coordinates, pivots, pivot_inverse)
end

DissipativeFrame(operators::Tuple) = build_dissipative_frame(operators)
DissipativeFrame(operators::AbstractVector) = build_dissipative_frame(operators)
DissipativeFrame(operators::SQA.QField...) = build_dissipative_frame(operators)

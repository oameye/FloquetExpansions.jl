# Dissipative-frame ownership and automatic discovery for positive completion.

Base.copy(frame::DissipativeFrame) = deepcopy(frame)

function Base.:(==)(
  left::DissipativeFrame{<:Tuple}, right::DissipativeFrame{<:Vector{SQA.QAdd}}
)
  return left.operators == Tuple(right.operators)
end
function Base.:(==)(
  left::DissipativeFrame{<:Vector{SQA.QAdd}}, right::DissipativeFrame{<:Tuple}
)
  return Tuple(left.operators) == right.operators
end
function Base.isequal(
  left::DissipativeFrame{<:Tuple}, right::DissipativeFrame{<:Vector{SQA.QAdd}}
)
  return isequal(left.operators, Tuple(right.operators))
end
function Base.isequal(
  left::DissipativeFrame{<:Vector{SQA.QAdd}}, right::DissipativeFrame{<:Tuple}
)
  return isequal(Tuple(left.operators), right.operators)
end
function Base.hash(frame::DissipativeFrame{<:Tuple}, h::UInt)
  return hash(:DissipativeFrame, hash(frame.operators, h))
end
function Base.hash(frame::DissipativeFrame{<:Vector{SQA.QAdd}}, h::UInt)
  return hash(:DissipativeFrame, hash(Tuple(frame.operators), h))
end

function coordinate_columns_independent(coordinates::KossakowskiMatrix)
  _, direction_count = size(coordinates)
  return length(coordinate_pivot_rows(coordinates)) == direction_count
end

function append_frame_candidate!(operators::Vector{SQA.QAdd}, operator::SQA.QField)
  iszero(dissipator(operator)) && return operators
  projected = projected_operator(operator)
  push!(operators, projected)
  coordinate_columns_independent(frame_direction_coordinates(operators)) || pop!(operators)
  return operators
end

function append_provenance_directions!(operators::Vector{SQA.QAdd}, ::NoProvenance)
  return operators
end

function append_provenance_directions!(
  operators::Vector{SQA.QAdd}, provenance::MicroscopicProvenance
)
  for seed in provenance.order
    operator = if seed.kind == COLLAPSE_SEED
      provenance.collapse_operators[seed.index]
    else
      provenance.jump_operators[seed.index]
    end
    append_frame_candidate!(operators, operator)
  end
  return operators
end

function append_generated_directions!(operators::Vector{SQA.QAdd}, L::Liouvillian)
  for operator in dissipative_support_directions(L)
    append_frame_candidate!(operators, operator)
  end
  return operators
end

function build_automatic_dissipative_frame(
  operators::Vector{SQA.QAdd}
)::DissipativeFrame{Vector{SQA.QAdd}}
  frame_operators = copy(operators)
  monomials, coordinates = frame_coordinate_data(frame_operators)
  pivots = independent_pivot_rows(coordinates)
  pivot_matrix = coefficient_matrix(length(pivots), length(frame_operators))
  for row in eachindex(pivots), column in eachindex(frame_operators)
    pivot_matrix[row, column] = coordinates[pivots[row], column]
  end
  pivot_inverse = inverse_coefficients(pivot_matrix)
  return DissipativeFrame(frame_operators, monomials, coordinates, pivots, pivot_inverse)
end

function automatic_dissipative_frame(
  expansion::FloquetExpansion
)::DissipativeFrame{Vector{SQA.QAdd}}
  operators = SQA.QAdd[]
  provenance = getfield(expansion, :provenance)
  append_provenance_directions!(operators, provenance)
  for component in getfield(expansion, :effective_components)
    append_generated_directions!(operators, component)
  end
  isempty(operators) && throw(
    ArgumentError(
      "positive completion found no dissipative directions; pass an explicit DissipativeFrame if a representation is required",
    ),
  )
  return build_automatic_dissipative_frame(operators)
end

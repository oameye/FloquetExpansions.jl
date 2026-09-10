# Dissipative-frame ownership and automatic discovery for positive completion.

Base.copy(frame::DissipativeFrame) = deepcopy(frame)

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

function automatic_dissipative_frame(expansion::FloquetExpansion)
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
  return DissipativeFrame(Tuple(operators))
end

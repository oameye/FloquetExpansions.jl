# Dissipative-frame ownership and automatic discovery for positive completion.

Base.copy(frame::DissipativeFrame) = deepcopy(frame)

function frame_direction_coordinates(operators::Vector{SQA.QAdd})
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

function coordinate_columns_independent(coordinates::KossakowskiMatrix)
  monomial_count, direction_count = size(coordinates)
  direction_count == 0 && return true
  monomial_count < direction_count && return false

  work = coefficient_matrix(direction_count, monomial_count)
  for row in 1:direction_count, column in 1:monomial_count
    work[row, column] = coordinates[column, row]
  end

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
      for trailing in column:monomial_count
        work[pivot_row, trailing], work[candidate, trailing] = work[candidate, trailing],
        work[pivot_row, trailing]
      end
    end

    pivot = work[pivot_row, column]
    for row in (pivot_row + 1):direction_count
      entry = simplify_coefficient(work[row, column])
      iszero(entry) && continue
      factor = simplify_coefficient(entry / pivot)
      for trailing in column:monomial_count
        work[row, trailing] = simplify_coefficient(
          work[row, trailing] - factor * work[pivot_row, trailing]
        )
      end
    end

    pivot_row += 1
    pivot_row > direction_count && return true
  end
  return false
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

function generated_support_directions(L::Liouvillian)
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

function raw_kossakowski_series(expansion::FloquetExpansion, frame::DissipativeFrame)
  components = getfield(expansion, :effective_components)
  result = KossakowskiMatrix[]
  sizehint!(result, length(components))
  for component in components
    push!(result, kossakowski(component, frame))
  end
  return result
end

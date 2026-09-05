function throw_recursive_elimination_obstruction(
  elimination::HermitianElimination, rate_order::Int
)
  if elimination.status == HERMITIAN_NEGATIVE_PIVOT
    throw(
      CompletionObstruction(
        rate_order,
        coefficient_from_completion(elimination.obstruction),
        :negative_direction,
      ),
    )
  elseif elimination.status == HERMITIAN_NONPOSITIVE_PIVOT
    throw(
      CompletionObstruction(
        rate_order,
        coefficient_from_completion(elimination.obstruction),
        :nonpositive_direction,
      ),
    )
  elseif elimination.status == HERMITIAN_ZERO_DIAGONAL_COUPLING
    throw(
      CompletionObstruction(
        rate_order,
        coefficient_from_completion(elimination.obstruction),
        :zero_diagonal_coupling,
      ),
    )
  end
  return elimination
end

function unresolved_recursive_leading_form(rate_order::Int)
  return CompletionObstruction(rate_order, convert(SQA.CNum, 0), :unresolved_leading_form)
end

function horizontal_series_stack(left::MatrixSeries, right::MatrixSeries, N::Int)
  validate_series_order(N)
  left_rows, left_columns = validate_matrix_series(left)
  right_rows, right_columns = validate_matrix_series(right)
  left_rows == right_rows ||
    throw(DimensionMismatch("series factors must have the same number of rows"))
  result = [completion_matrix_zeros(left_rows, left_columns + right_columns) for _ in 0:N]
  for grade in 0:N
    result[grade + 1] = hcat(
      matrix_coefficient(left, grade, left_rows, left_columns),
      matrix_coefficient(right, grade, right_rows, right_columns),
    )
  end
  return result
end

function upper_zero_pad_series(factor::MatrixSeries, rows::Int, N::Int)
  validate_series_order(N)
  factor_rows, columns = validate_matrix_series(factor)
  upper = [completion_matrix_zeros(rows, columns) for _ in 0:N]
  lower = [matrix_coefficient(factor, grade, factor_rows, columns) for grade in 0:N]
  return vertical_series_stack(upper, lower, N)
end

function classify_recursive_onset!(
  series::MatrixSeries, onset::Int, conditions::CompletionConditions, rate_offset::Int
)
  leading = hermitian_eliminate(series[onset + 1], conditions)
  global_rate_order = rate_offset + onset
  throw_recursive_elimination_obstruction(leading, global_rate_order)
  leading.active_rank > 0 || throw(unresolved_recursive_leading_form(global_rate_order))
  isodd(onset) && throw(FractionalJumpOnset(global_rate_order))
  return leading
end

function recurse_from_onset(
  series::MatrixSeries,
  onset::Int,
  N::Int,
  conditions::CompletionConditions,
  rate_offset::Int,
)
  classify_recursive_onset!(series, onset, conditions, rate_offset)
  residual_order = N - onset
  shifted = series_shift_down(series, onset, residual_order)
  factor, stages = recursive_gram_factor_series(
    shifted, residual_order, conditions, rate_offset + onset
  )
  return series_scale_shift(factor, onset ÷ 2, N), stages
end

function recursive_gram_factor_series(
  series::MatrixSeries, N::Int, conditions::CompletionConditions, rate_offset::Int=0
)
  q, columns = validate_matrix_series(series)
  q == columns || throw(DimensionMismatch("Kossakowski series must be square"))
  hermitian_series(series) ||
    throw(ArgumentError("Gram completion requires a Hermitian retained Kossakowski series"))

  onset = matrix_series_onset(series, N)
  if onset < 0
    factor = [completion_matrix_zeros(q, 0) for _ in 0:N]
    return factor, GramStage[GramStage(rate_offset, 0, q)]
  elseif onset > 0
    factor, stages = recurse_from_onset(series, onset, N, conditions, rate_offset)
    return factor, vcat(GramStage[GramStage(rate_offset, 0, q)], stages)
  end

  elimination = hermitian_eliminate(series[1], conditions)
  throw_recursive_elimination_obstruction(elimination, rate_offset)
  active_rank = elimination.active_rank
  dark_rank = q - active_rank
  active_rank > 0 || throw(unresolved_recursive_leading_form(rate_offset))

  reduced = apply_congruence(series, elimination.transform, N)
  active = matrix_series_block(reduced, 1:active_rank, 1:active_rank)
  ldl = graded_ldl(active, N, conditions)
  active_factor = ldl_gram_factor(ldl, N, conditions)
  stages = GramStage[GramStage(rate_offset, active_rank, dark_rank)]

  reduced_factor = if dark_rank == 0
    active_factor
  else
    dark_columns = (active_rank + 1):q
    cross = matrix_series_block(reduced, 1:active_rank, dark_columns)
    dark = matrix_series_block(reduced, dark_columns, dark_columns)
    residual = series_schur(active, cross, dark, N, conditions)

    solved = series_solve(active_factor, cross, N, conditions)
    dark_rows = series_adjoint(solved)
    dressed_active = vertical_series_stack(active_factor, dark_rows, N)

    residual_onset = matrix_series_onset(residual, N)
    if residual_onset < 0
      dressed_active
    else
      residual_factor, residual_stages = recurse_from_onset(
        residual, residual_onset, N, conditions, rate_offset
      )
      append!(stages, residual_stages)
      embedded_residual = upper_zero_pad_series(residual_factor, active_rank, N)
      horizontal_series_stack(dressed_active, embedded_residual, N)
    end
  end

  factor = undo_congruence_factor(reduced_factor, elimination.transform, N, conditions)
  return factor, stages
end

function recursive_gram_positive_completion(
  expansion::FloquetExpansion, frame::DissipativeFrame, algorithm::Gram
)
  N = getfield(expansion, :order) - 1
  raw_matrices = raw_kossakowski_series(expansion, frame)
  series = completion_series(raw_matrices)
  conditions = CompletionConditions()
  seed_completion_conditions!(conditions, getfield(expansion, :provenance))

  factor, stages = recursive_gram_factor_series(series, N, conditions)
  amplitudes = physical_factor_amplitudes(factor, getfield(expansion, :generator).wd, N)
  finite = finite_factor(amplitudes)
  completed_matrix = multiply_coefficients(finite, adjoint_coefficients(finite))
  simplify_matrix!(completed_matrix)
  matrix_is_hermitian(completed_matrix) ||
    throw(ArgumentError("completed Gram Kossakowski matrix is not Hermitian"))

  completed_channels = completed_collapse_channels(frame, finite)
  raw_generator = effective_generator(expansion)
  coherent = hamiltonian(raw_generator, frame)
  generator = liouvillian(coherent; channels=completed_channels)

  factorization = GramFactorization(amplitudes, factor_onsets(amplitudes), stages)
  completion = PositiveCompletion(
    algorithm,
    frame,
    completed_matrix,
    completed_channels,
    condition_coefficients(conditions.positivity),
    condition_coefficients(conditions.regularity),
    factorization,
    generator,
  )
  return with_completion(expansion, completion)
end

function positive_completion_impl(
  expansion::FloquetExpansion{G,P,E,Uncompleted,R}, algorithm::Gram
) where {G,P,E,R}
  return recursive_gram_positive_completion(
    expansion, automatic_dissipative_frame(expansion), algorithm
  )
end

function positive_completion_impl(
  expansion::FloquetExpansion{G,P,E,Uncompleted,R}, algorithm::Gram, frame::DissipativeFrame
) where {G,P,E,R}
  return recursive_gram_positive_completion(expansion, frame, algorithm)
end

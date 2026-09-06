"""
    GramStage

Compact diagnostic record for one algebraic Gram-completion stage. `grade` is the
inverse-frequency grade at which the stage is opened, while `active_rank` and `dark_rank`
record the active and dark dimensions of the Hermitian leading form.
"""
struct GramStage
  grade::Int
  active_rank::Int
  dark_rank::Int
end

"""
    GramFactorization <: CompletionFactorization

Factorization data produced by [`Gram`](@ref) positive completion.

`amplitudes[n + 1]` is the physical order-`n` collapse-amplitude matrix, including its
`ωd^(-n)` scaling. All amplitude matrices have the same number of columns. `onsets` gives
the first nonzero grade of each completed collapse channel and `stages` records the compact
active/dark elimination history.
"""
struct GramFactorization <: CompletionFactorization
  amplitudes::Vector{KossakowskiMatrix}
  onsets::Vector{Int}
  stages::Vector{GramStage}
end

function Base.copy(factorization::GramFactorization)
  return GramFactorization(
    [copy(amplitude) for amplitude in factorization.amplitudes],
    copy(factorization.onsets),
    copy(factorization.stages),
  )
end

function matrix_series_block(
  series::MatrixSeries, rows::UnitRange{Int}, columns::UnitRange{Int}
)
  return [Matrix{CompletionScalar}(matrix[rows, columns]) for matrix in series]
end

function matrix_series_onset(series::MatrixSeries, N::Int)
  validate_series_order(N)
  for grade in 0:min(N, length(series) - 1)
    structurally_zero(series[grade + 1]) || return grade
  end
  return -1
end

function throw_elimination_obstruction(elimination::HermitianElimination, grade::Int)
  if elimination.status == HERMITIAN_NEGATIVE_PIVOT ||
    elimination.status == HERMITIAN_NONPOSITIVE_PIVOT
    throw(
      ArgumentError(
        "Gram completion encountered a retained negative Hermitian direction at grade $grade: $(elimination.obstruction)",
      ),
    )
  elseif elimination.status == HERMITIAN_ZERO_DIAGONAL_COUPLING
    throw(
      ArgumentError(
        "Gram completion encountered a zero-diagonal/nonzero-coupling PSD obstruction at grade $grade: $(elimination.obstruction)",
      ),
    )
  end
  return elimination
end

function classify_unresolved_dark_series!(
  residual::MatrixSeries, N::Int, conditions::CompletionConditions
)
  onset = matrix_series_onset(residual, N)
  onset < 0 && return nothing

  leading = hermitian_eliminate_structured(residual[onset + 1], conditions)
  throw_elimination_obstruction(leading, onset)
  leading.active_rank == 0 && throw(
    ArgumentError(
      "Gram completion could not resolve the dark-sector leading form at grade $onset"
    ),
  )

  isodd(onset) && throw(
    ArgumentError(
      "Gram completion requires a fractional/Puiseux collapse-amplitude onset at grade $onset; recursive dark-sector handling is required",
    ),
  )
  return throw(
    ArgumentError(
      "Gram completion found a new even dark-sector onset at grade $onset; recursive dark-sector handling is implemented in the next completion stage",
    ),
  )
end

function vertical_series_stack(upper::MatrixSeries, lower::MatrixSeries, N::Int)
  validate_series_order(N)
  upper_rows, columns = validate_matrix_series(upper)
  lower_rows, lower_columns = validate_matrix_series(lower)
  columns == lower_columns ||
    throw(DimensionMismatch("series factors must have the same number of columns"))
  result = [completion_matrix_zeros(upper_rows + lower_rows, columns) for _ in 0:N]
  for grade in 0:N
    result[grade + 1] = vcat(upper[grade + 1], lower[grade + 1])
  end
  return result
end

function gram_factor_series(series::MatrixSeries, N::Int, conditions::CompletionConditions)
  q, columns = validate_matrix_series(series)
  q == columns || throw(DimensionMismatch("Kossakowski series must be square"))
  hermitian_series(series) ||
    throw(ArgumentError("Gram completion requires a Hermitian retained Kossakowski series"))

  elimination = hermitian_eliminate_structured(series[1], conditions)
  throw_elimination_obstruction(elimination, 0)
  active_rank = elimination.active_rank
  dark_rank = q - active_rank

  if active_rank == 0
    matrix_series_onset(series, N) < 0 ||
      classify_unresolved_dark_series!(series, N, conditions)
    factor = [completion_matrix_zeros(q, 0) for _ in 0:N]
    return factor, GramStage(0, 0, q)
  end

  reduced = apply_congruence(series, elimination.transform, N)
  active = matrix_series_block(reduced, 1:active_rank, 1:active_rank)
  ldl = graded_ldl(active, N, conditions)
  active_factor = ldl_gram_factor(ldl, N, conditions)

  reduced_factor = if dark_rank == 0
    active_factor
  else
    dark_columns = (active_rank + 1):q
    cross = matrix_series_block(reduced, 1:active_rank, dark_columns)
    dark = matrix_series_block(reduced, dark_columns, dark_columns)
    residual, dark_rows = gram_feshbach_dressing(active_factor, cross, dark, N, conditions)
    classify_unresolved_dark_series!(residual, N, conditions)
    vertical_series_stack(active_factor, dark_rows, N)
  end

  factor = planned_undo_congruence_factor(
    reduced_factor, elimination.transform, N, conditions
  )
  return factor, GramStage(0, active_rank, dark_rank)
end

function physical_factor_amplitudes(factor::MatrixSeries, wd::Symbolics.Num, N::Int)
  rows, columns = validate_matrix_series(factor)
  amplitudes = KossakowskiMatrix[]
  sizehint!(amplitudes, N + 1)
  for grade in 0:N
    scale = completion_scalar(iszero(grade) ? 1 : wd^(-grade))
    matrix = completion_matrix_zeros(rows, columns)
    for index in eachindex(matrix)
      matrix[index] = simplify_scalar(scale * factor[grade + 1][index])
    end
    push!(amplitudes, coefficient_matrix_from_completion(matrix))
  end
  return amplitudes
end

function factor_onsets(amplitudes::Vector{KossakowskiMatrix})
  isempty(amplitudes) && return Int[]
  channels = size(first(amplitudes), 2)
  result = fill(-1, channels)
  for channel in 1:channels
    for grade in eachindex(amplitudes)
      any(!iszero(amplitudes[grade][row, channel]) for row in axes(amplitudes[grade], 1)) ||
        continue
      result[channel] = grade - 1
      break
    end
  end
  return result
end

function finite_factor(amplitudes::Vector{KossakowskiMatrix})
  isempty(amplitudes) && return coefficient_matrix(0, 0)
  rows, columns = size(first(amplitudes))
  result = coefficient_matrix(rows, columns)
  for amplitude in amplitudes
    size(amplitude) == (rows, columns) ||
      throw(DimensionMismatch("all Gram amplitudes must have the same dimensions"))
    for index in eachindex(result)
      result[index] = simplify_coefficient(result[index] + amplitude[index])
    end
  end
  return result
end

function completed_collapse_channels(frame::DissipativeFrame, factor::KossakowskiMatrix)
  size(factor, 1) == length(frame.operators) ||
    throw(DimensionMismatch("Gram factor does not match the dissipative frame"))
  result = CollapseChannel{SQA.QAdd}[]
  for column in axes(factor, 2)
    operator = zero(SQA.QAdd)
    for row in axes(factor, 1)
      coefficient = factor[row, column]
      iszero(coefficient) && continue
      operator = operator + coefficient * frame.operators[row]
    end
    operator = SQA.simplify(operator)
    iszero(operator) || push!(result, collapse(operator))
  end
  return result
end

function gram_positive_completion(
  expansion::FloquetExpansion, frame::DissipativeFrame, algorithm::Gram
)
  N = getfield(expansion, :order) - 1
  raw_matrices = raw_kossakowski_series(expansion, frame)
  series = completion_series(raw_matrices)
  conditions = CompletionConditions()
  seed_completion_conditions!(conditions, getfield(expansion, :provenance))

  factor, stage = gram_factor_series(series, N, conditions)
  amplitudes = physical_factor_amplitudes(factor, getfield(expansion, :generator).wd, N)
  finite = finite_factor(amplitudes)
  completed_matrix = multiply_coefficients(finite, adjoint_coefficients(finite))
  simplify_matrix!(completed_matrix)
  completed_channels = completed_collapse_channels(frame, finite)
  factorization = GramFactorization(amplitudes, factor_onsets(amplitudes), [stage])

  return finalize_positive_completion(
    expansion,
    algorithm,
    frame,
    raw_matrices,
    completed_matrix,
    completed_channels,
    conditions,
    factorization,
  )
end

"""
    dissipative_frame(expansion::FloquetExpansion)

Return an independent copy of the ordered dissipative frame owned by a positively completed
Floquet expansion.
"""
function dissipative_frame(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:PositiveCompletion,R}
  return copy(stored_dissipative_frame(expansion))
end

"""
    channels(expansion::FloquetExpansion)

Return the completed physical channels. Together with [`hamiltonian`](@ref), these channels
reconstruct [`effective_generator`](@ref) exactly.
"""
function channels(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:PositiveCompletion,R}
  return deepcopy(stored_completion(expansion).channels)
end

"""
    positivity_conditions(expansion::FloquetExpansion)

Return scalar conditions assumed nonnegative on the symbolic stratum used by positive
completion.
"""
function positivity_conditions(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:PositiveCompletion,R}
  return copy(stored_completion(expansion).positivity_conditions)
end

"""
    regularity_conditions(expansion::FloquetExpansion)

Return scalar conditions assumed nonzero to remain on the fixed-rank symbolic stratum used by
positive completion.
"""
function regularity_conditions(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:PositiveCompletion,R}
  return copy(stored_completion(expansion).regularity_conditions)
end

"""
    factorization(expansion::FloquetExpansion)

Return the algorithm-specific factorization data stored by a positively completed expansion.
For [`Gram`](@ref) completion this is a [`GramFactorization`](@ref).
"""
function factorization(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:PositiveCompletion,R}
  return copy(stored_completion(expansion).factorization)
end

"""
    kossakowski(expansion::FloquetExpansion)

Return the finite completed Kossakowski matrix in [`dissipative_frame`](@ref). This no-frame
form is defined for positively completed Liouvillian expansions.
"""
function kossakowski(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian,C<:PositiveCompletion,R}
  return copy(stored_completion(expansion).kossakowski)
end

"""
    kossakowski_component(expansion::FloquetExpansion, n::Int)

Return the cached retained order-`n` Kossakowski contribution in the completed expansion's
stored dissipative frame. Positive completion does not alter retained components.
"""
function kossakowski_component(
  expansion::FloquetExpansion{G,P,E,C,R}, n::Int
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian,C<:PositiveCompletion,R}
  retained = stored_retained_kossakowski(expansion)
  0 <= n < length(retained) || throw(BoundsError(retained, n + 1))
  return copy(retained[n + 1])
end

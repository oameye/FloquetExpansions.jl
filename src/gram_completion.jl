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

completion_scalar(x::SQA.CNum) = completion_scalar(SQA.to_num(x))

function completion_matrix(matrix::KossakowskiMatrix)
  result = completion_matrix_zeros(size(matrix, 1), size(matrix, 2))
  for index in eachindex(matrix)
    result[index] = completion_scalar(matrix[index])
  end
  return result
end

function coefficient_from_completion(value::CompletionScalar)
  return SQA.simplify(convert(SQA.CNum, simplify_scalar(value)))
end

function coefficient_matrix_from_completion(matrix::CompletionMatrix)
  result = coefficient_matrix(size(matrix, 1), size(matrix, 2))
  for index in eachindex(matrix)
    result[index] = coefficient_from_completion(matrix[index])
  end
  return result
end

function condition_coefficients(values::Vector{CompletionScalar})
  result = SQA.CNum[]
  sizehint!(result, length(values))
  for value in values
    push!(result, coefficient_from_completion(value))
  end
  return result
end

function seed_completion_conditions!(conditions::CompletionConditions, ::NoProvenance)
  return conditions
end

function seed_completion_conditions!(
  conditions::CompletionConditions, provenance::MicroscopicProvenance
)
  for assumption in provenance.rate_assumptions
    require_positivity!(conditions, completion_scalar(assumption.rate))
  end
  return conditions
end

function dependent_frame_error(error::ArgumentError)
  return occursin("linearly dependent modulo the identity", sprint(showerror, error))
end

function append_frame_candidate!(operators::Vector{SQA.QAdd}, operator::SQA.QField)
  iszero(dissipator(operator)) && return operators
  projected = projected_operator(operator)
  if isempty(operators)
    push!(operators, projected)
    return operators
  end

  trial = (Tuple(operators)..., projected)
  try
    DissipativeFrame(trial)
  catch error
    error isa ArgumentError || rethrow()
    dependent_frame_error(error) || rethrow()
    return operators
  end
  push!(operators, projected)
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
  has_two_sided_terms(L) || return operators
  generated = support_frame(L)
  for operator in generated.operators
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

function completion_series(matrices::Vector{KossakowskiMatrix})
  result = CompletionMatrix[]
  sizehint!(result, length(matrices))
  for matrix in matrices
    push!(result, completion_matrix(matrix))
  end
  return result
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

  leading = hermitian_eliminate(residual[onset + 1], conditions)
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

  elimination = hermitian_eliminate(series[1], conditions)
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
    residual = series_schur(active, cross, dark, N, conditions)
    classify_unresolved_dark_series!(residual, N, conditions)

    solved = series_solve(active_factor, cross, N, conditions)
    dark_rows = series_adjoint(solved)
    vertical_series_stack(active_factor, dark_rows, N)
  end

  factor = undo_congruence_factor(reduced_factor, elimination.transform, N, conditions)
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
  matrix_is_hermitian(completed_matrix) ||
    throw(ArgumentError("completed Gram Kossakowski matrix is not Hermitian"))

  completed_channels = completed_collapse_channels(frame, finite)
  raw_generator = effective_generator(expansion)
  coherent = hamiltonian(raw_generator, frame)
  generator = liouvillian(coherent; channels=completed_channels)

  factorization = GramFactorization(amplitudes, factor_onsets(amplitudes), [stage])
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

function positive_completion_impl(expansion::FloquetExpansion, algorithm::Gram)
  return gram_positive_completion(
    expansion, automatic_dissipative_frame(expansion), algorithm
  )
end

function positive_completion_impl(
  expansion::FloquetExpansion, algorithm::Gram, frame::DissipativeFrame
)
  return gram_positive_completion(expansion, frame, algorithm)
end

"""
    dissipative_frame(expansion::FloquetExpansion)

Return the ordered dissipative frame stored by a positively completed Floquet expansion.
"""
function dissipative_frame(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:PositiveCompletion,R}
  return getfield(expansion, :completion).frame
end

"""
    channels(expansion::FloquetExpansion)

Return the completed physical collapse channels. Together with [`hamiltonian`](@ref), these
channels reconstruct [`effective_generator`](@ref) exactly.
"""
function channels(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:PositiveCompletion,R}
  return getfield(expansion, :completion).channels
end

"""
    positivity_conditions(expansion::FloquetExpansion)

Return scalar conditions assumed nonnegative on the symbolic stratum used by positive
completion.
"""
function positivity_conditions(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:PositiveCompletion,R}
  return copy(getfield(expansion, :completion).positivity_conditions)
end

"""
    regularity_conditions(expansion::FloquetExpansion)

Return scalar conditions assumed nonzero to remain on the fixed-rank symbolic stratum used by
positive completion.
"""
function regularity_conditions(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:PositiveCompletion,R}
  return copy(getfield(expansion, :completion).regularity_conditions)
end

"""
    factorization(expansion::FloquetExpansion)

Return the algorithm-specific factorization data stored by a positively completed expansion.
For [`Gram`](@ref) completion this is a [`GramFactorization`](@ref).
"""
function factorization(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:PositiveCompletion,R}
  return getfield(expansion, :completion).factorization
end

"""
    kossakowski(expansion::FloquetExpansion)

Return the finite completed Kossakowski matrix in [`dissipative_frame`](@ref). This no-frame
form is defined for positively completed Liouvillian expansions.
"""
function kossakowski(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian,C<:PositiveCompletion,R}
  return getfield(expansion, :completion).kossakowski
end

"""
    kossakowski_component(expansion::FloquetExpansion, n::Int)

Return the retained order-`n` Kossakowski contribution in the completed expansion's stored
[`dissipative_frame`](@ref). Positive completion does not alter retained components.
"""
function kossakowski_component(
  expansion::FloquetExpansion{G,P,E,C,R}, n::Int
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian,C<:PositiveCompletion,R}
  return kossakowski_component(expansion, dissipative_frame(expansion), n)
end

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
  return copy(stored_completion(expansion).channels)
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

"""
    SpectralFactorization <: CompletionFactorization

Diagnostic data produced by [`Spectral`](@ref) positive completion. `rates` are the finite
HCM-completed branch rates, `vectors` are the corresponding normalized dissipative-frame
coordinates, `onsets` record the first retained rate order of each branch, and `puiseux`
marks odd rate onsets whose rate-folded collapse amplitudes would begin at half-integer
order.
"""
struct SpectralFactorization <: CompletionFactorization
  rates::Vector{SQA.CNum}
  vectors::Vector{Vector{SQA.CNum}}
  onsets::Vector{Int}
  puiseux::Vector{Bool}
end

function spectral_leading_diagonal(series::MatrixSeries)
  q, columns = validate_matrix_series(series)
  q == columns || throw(DimensionMismatch("Kossakowski series must be square"))
  hermitian_series(series) ||
    throw(ArgumentError("Spectral completion requires a Hermitian retained Kossakowski series"))

  leading = series[1]
  for column in 1:q, row in 1:q
    row == column && continue
    structurally_zero(leading[row, column]) || throw(
      ArgumentError(
        "Spectral completion requires d^(0) to be diagonal in the supplied dissipative frame",
      ),
    )
  end
  return [hermitian_real(leading[index, index]) for index in 1:q]
end

function completion_vector_zeros(length::Int)
  length >= 0 || throw(ArgumentError("vector length must be nonnegative"))
  return fill(completion_zero(), length)
end

function completion_matvec(matrix::CompletionMatrix, vector::Vector{CompletionScalar})
  size(matrix, 2) == length(vector) || throw(DimensionMismatch("matrix/vector mismatch"))
  result = completion_vector_zeros(size(matrix, 1))
  for row in axes(matrix, 1)
    value = completion_zero()
    for column in axes(matrix, 2)
      value += matrix[row, column] * vector[column]
    end
    result[row] = simplify_scalar(value)
  end
  return result
end

function completion_dot(left::Vector{CompletionScalar}, right::Vector{CompletionScalar})
  length(left) == length(right) || throw(DimensionMismatch("vector dimensions must match"))
  result = completion_zero()
  for index in eachindex(left)
    result += conj(left[index]) * right[index]
  end
  return simplify_scalar(result)
end

function spectral_branch_series(
  series::MatrixSeries,
  leading_rates::Vector{CompletionScalar},
  branch::Int,
  N::Int,
  conditions::CompletionConditions,
)
  q = length(leading_rates)
  vectors = [completion_vector_zeros(q) for _ in 0:N]
  vectors[1][branch] = completion_one()
  rates = fill(completion_zero(), N + 1)
  rates[1] = leading_rates[branch]

  for order in 1:N
    rate = completion_zero()
    for k in 1:order
      k + 1 <= length(series) || continue
      contribution = completion_matvec(series[k + 1], vectors[order - k + 1])
      rate += contribution[branch]
    end
    for k in 1:(order - 1)
      rate -= rates[k + 1] * vectors[order - k + 1][branch]
    end
    rates[order + 1] = hermitian_real(simplify_scalar(rate))

    residual = completion_vector_zeros(q)
    for k in 1:order
      k + 1 <= length(series) || continue
      contribution = completion_matvec(series[k + 1], vectors[order - k + 1])
      for index in 1:q
        residual[index] += contribution[index]
      end
    end
    for k in 1:order
      for index in 1:q
        residual[index] -= rates[k + 1] * vectors[order - k + 1][index]
      end
    end

    for index in 1:q
      index == branch && continue
      value = simplify_scalar(residual[index])
      gap = simplify_scalar(leading_rates[branch] - leading_rates[index])
      if structurally_zero(gap)
        structurally_zero(value) || throw(
          ArgumentError(
            "Spectral completion encountered unresolved mixing inside a degenerate leading sector at order $order",
          ),
        )
        continue
      end
      structurally_nonzero(gap, conditions) || require_regularity!(conditions, gap)
      vectors[order + 1][index] = simplify_scalar(value / gap)
    end

    normalization = completion_zero()
    for k in 1:(order - 1)
      normalization += completion_dot(vectors[k + 1], vectors[order - k + 1])
    end
    vectors[order + 1][branch] =
      simplify_scalar(-hermitian_real(normalization) / completion_scalar(2))
  end

  return rates, vectors
end

function hcm_completed_rate(
  rates::ScalarSeries,
  wd::Symbolics.Num,
  N::Int,
  conditions::CompletionConditions,
)
  onset = series_onset(rates, N)
  onset < 0 && return coefficient_zero(), -1, false

  leading = hermitian_real(rates[onset + 1])
  sign = structural_sign(leading, conditions)
  if sign == SIGN_NEGATIVE || sign == SIGN_NONPOSITIVE
    throw(
      CompletionObstruction(onset, coefficient_from_completion(leading), :negative_spectral_rate)
    )
  elseif sign == SIGN_UNKNOWN
    require_positivity!(conditions, leading)
  end

  retained = N - onset
  shifted = series_shift_down(rates, onset, retained)
  root = scalar_series_sqrt(shifted, retained, conditions)
  finite_root = completion_zero()
  for grade in 0:retained
    scale = completion_scalar(iszero(grade) ? 1 : wd^(-grade))
    finite_root += scale * root[grade + 1]
  end
  rate = completion_scalar(iszero(onset) ? 1 : wd^(-onset)) * finite_root^2
  return coefficient_from_completion(hermitian_real(simplify_scalar(rate))), onset, isodd(onset)
end

function finite_spectral_vector(
  vectors::Vector{Vector{CompletionScalar}}, wd::Symbolics.Num
)
  q = length(first(vectors))
  finite = completion_vector_zeros(q)
  for grade in 0:(length(vectors) - 1)
    scale = completion_scalar(iszero(grade) ? 1 : wd^(-grade))
    for index in 1:q
      finite[index] += scale * vectors[grade + 1][index]
    end
  end

  norm_squared = hermitian_real(completion_dot(finite, finite))
  structurally_zero(norm_squared) &&
    throw(ArgumentError("Spectral completion produced a zero perturbative branch vector"))
  norm = completion_scalar(sqrt(real(norm_squared)))
  return [coefficient_from_completion(simplify_scalar(value / norm)) for value in finite]
end

function spectral_completed_matrix(
  rates::Vector{SQA.CNum}, vectors::Vector{Vector{SQA.CNum}}, q::Int
)
  result = coefficient_matrix(q, q)
  for branch in eachindex(rates)
    iszero(rates[branch]) && continue
    vector = vectors[branch]
    for row in 1:q, column in 1:q
      result[row, column] = simplify_coefficient(
        result[row, column] + rates[branch] * vector[row] * conj(vector[column])
      )
    end
  end
  return simplify_matrix!(result)
end

function spectral_completed_channels(
  frame::DissipativeFrame,
  rates::Vector{SQA.CNum},
  vectors::Vector{Vector{SQA.CNum}},
)
  result = RateWeightedJump{SQA.QAdd}[]
  for branch in eachindex(rates)
    rate = rates[branch]
    iszero(rate) && continue
    operator = zero(SQA.QAdd)
    for index in eachindex(frame.operators)
      coefficient = vectors[branch][index]
      iszero(coefficient) && continue
      operator = operator + coefficient * frame.operators[index]
    end
    operator = SQA.simplify(operator)
    iszero(operator) || push!(result, jump(operator, rate))
  end
  return result
end

function spectral_positive_completion(
  expansion::FloquetExpansion, frame::DissipativeFrame, algorithm::Spectral
)
  N = getfield(expansion, :order) - 1
  raw_matrices = raw_kossakowski_series(expansion, frame)
  series = completion_series(raw_matrices)
  conditions = CompletionConditions()
  seed_completion_conditions!(conditions, getfield(expansion, :provenance))
  leading_rates = spectral_leading_diagonal(series)

  wd = getfield(expansion, :generator).wd
  completed_rates = SQA.CNum[]
  completed_vectors = Vector{SQA.CNum}[]
  onsets = Int[]
  puiseux = Bool[]
  for branch in eachindex(leading_rates)
    rates, vectors = spectral_branch_series(series, leading_rates, branch, N, conditions)
    completed_rate, onset, fractional = hcm_completed_rate(rates, wd, N, conditions)
    push!(completed_rates, completed_rate)
    push!(completed_vectors, finite_spectral_vector(vectors, wd))
    push!(onsets, onset)
    push!(puiseux, fractional)
  end

  completed_matrix = spectral_completed_matrix(
    completed_rates, completed_vectors, length(frame.operators)
  )
  matrix_is_hermitian(completed_matrix) ||
    throw(ArgumentError("completed spectral Kossakowski matrix is not Hermitian"))
  completed_channels = spectral_completed_channels(frame, completed_rates, completed_vectors)

  raw_generator = effective_generator(expansion)
  coherent = hamiltonian(raw_generator, frame)
  generator = liouvillian(coherent; channels=completed_channels)
  factorization = SpectralFactorization(completed_rates, completed_vectors, onsets, puiseux)
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

function positive_completion_impl(expansion::FloquetExpansion, algorithm::Spectral)
  return spectral_positive_completion(
    expansion, automatic_dissipative_frame(expansion), algorithm
  )
end

function positive_completion_impl(
  expansion::FloquetExpansion, algorithm::Spectral, frame::DissipativeFrame
)
  return spectral_positive_completion(expansion, frame, algorithm)
end

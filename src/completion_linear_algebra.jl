# Structured symbolic linear algebra used by the positive-completion kernels.
#
# Keep this separate from the coefficient/series recurrences in `matrix_series.jl`: the
# routines below cache structural work that depends only on a leading matrix and exploit
# Hermitian/triangular structure without materializing dense elementary transforms.

struct MatrixSolvePlan
  lower::CompletionMatrix
  upper::CompletionMatrix
  permutation::Vector{Int}
  inverse_diagonal::Vector{CompletionScalar}
end

function matrix_solve_plan(A::CompletionMatrix, conditions::CompletionConditions)
  n, m = size(A)
  n == m || throw(DimensionMismatch("solve matrix must be square"))

  lower = completion_identity(n)
  upper = copy(A)
  permutation = collect(1:n)

  for k in 1:n
    pivot_row = choose_pivot_row(upper, k, conditions)
    pivot_row == 0 &&
      throw(ArgumentError("matrix is singular on the current symbolic stratum"))

    if pivot_row != k
      for column in 1:n
        upper[k, column], upper[pivot_row, column] = upper[pivot_row, column],
        upper[k, column]
      end
      for column in 1:(k - 1)
        lower[k, column], lower[pivot_row, column] = lower[pivot_row, column],
        lower[k, column]
      end
      permutation[k], permutation[pivot_row] = permutation[pivot_row], permutation[k]
    end

    pivot = simplify_scalar(upper[k, k])
    structurally_nonzero(pivot, conditions) || require_regularity!(conditions, pivot)
    upper[k, k] = pivot

    for row in (k + 1):n
      entry = simplify_scalar(upper[row, k])
      structurally_zero(entry) && continue
      factor = simplify_scalar(entry / pivot)
      lower[row, k] = factor
      upper[row, k] = completion_zero()
      for column in (k + 1):n
        upper[row, column] = simplify_scalar(upper[row, column] - factor * upper[k, column])
      end
    end
  end

  inverse_diagonal = Vector{CompletionScalar}(undef, n)
  for index in 1:n
    pivot = simplify_scalar(upper[index, index])
    structurally_nonzero(pivot, conditions) || require_regularity!(conditions, pivot)
    inverse_diagonal[index] = simplify_scalar(completion_one() / pivot)
  end

  return MatrixSolvePlan(lower, upper, permutation, inverse_diagonal)
end

function apply_solve_plan(plan::MatrixSolvePlan, B::CompletionMatrix)
  n = size(plan.upper, 1)
  size(B, 1) == n || throw(DimensionMismatch("right-hand side has incompatible row count"))

  rhs = completion_matrix_zeros(n, size(B, 2))
  for column in axes(B, 2), row in 1:n
    rhs[row, column] = B[plan.permutation[row], column]
  end
  forward = unit_lower_triangular_solve(plan.lower, rhs)

  result = completion_matrix_zeros(n, size(B, 2))
  for column in axes(B, 2), row in n:-1:1
    residual = forward[row, column]
    for trailing in (row + 1):n
      residual -= plan.upper[row, trailing] * result[trailing, column]
    end
    result[row, column] = simplify_scalar(residual * plan.inverse_diagonal[row])
  end
  return result
end

function planned_series_solve(
  A::MatrixSeries, B::MatrixSeries, N::Int, conditions::CompletionConditions
)
  validate_series_order(N)
  n, m = validate_matrix_series(A)
  n == m || throw(DimensionMismatch("left matrix series must be square"))
  b_rows, b_columns = validate_matrix_series(B)
  b_rows == n ||
    throw(DimensionMismatch("matrix-series right-hand side has incompatible rows"))

  plan = matrix_solve_plan(A[1], conditions)
  result = [completion_matrix_zeros(n, b_columns) for _ in 0:N]
  for order in 0:N
    rhs = copy(matrix_coefficient(B, order, b_rows, b_columns))
    for k in 1:order
      k + 1 <= length(A) || continue
      rhs -= A[k + 1] * result[order - k + 1]
    end
    result[order + 1] = apply_solve_plan(plan, rhs)
  end
  return result
end

struct LowerTriangularSolvePlan
  lower::CompletionMatrix
  inverse_diagonal::Vector{CompletionScalar}
end

function lower_triangular_solve_plan(L::CompletionMatrix, conditions::CompletionConditions)
  n, m = size(L)
  n == m || throw(DimensionMismatch("triangular solve matrix must be square"))
  for row in 1:n, column in (row + 1):n
    structurally_zero(L[row, column]) ||
      throw(ArgumentError("matrix is not structurally lower triangular"))
  end

  inverse_diagonal = Vector{CompletionScalar}(undef, n)
  for index in 1:n
    pivot = simplify_scalar(L[index, index])
    structurally_nonzero(pivot, conditions) || require_regularity!(conditions, pivot)
    inverse_diagonal[index] = simplify_scalar(completion_one() / pivot)
  end
  return LowerTriangularSolvePlan(copy(L), inverse_diagonal)
end

function apply_solve_plan(plan::LowerTriangularSolvePlan, B::CompletionMatrix)
  n = size(plan.lower, 1)
  size(B, 1) == n || throw(DimensionMismatch("right-hand side has incompatible row count"))
  result = completion_matrix_zeros(n, size(B, 2))
  for column in axes(B, 2), row in 1:n
    residual = B[row, column]
    for previous in 1:(row - 1)
      residual -= plan.lower[row, previous] * result[previous, column]
    end
    result[row, column] = simplify_scalar(residual * plan.inverse_diagonal[row])
  end
  return result
end

function lower_triangular_series_solve(
  L::MatrixSeries, B::MatrixSeries, N::Int, conditions::CompletionConditions
)
  validate_series_order(N)
  n, m = validate_matrix_series(L)
  n == m || throw(DimensionMismatch("left matrix series must be square"))
  b_rows, b_columns = validate_matrix_series(B)
  b_rows == n ||
    throw(DimensionMismatch("matrix-series right-hand side has incompatible rows"))

  plan = lower_triangular_solve_plan(L[1], conditions)
  result = [completion_matrix_zeros(n, b_columns) for _ in 0:N]
  for order in 0:N
    rhs = copy(matrix_coefficient(B, order, b_rows, b_columns))
    for k in 1:order
      k + 1 <= length(L) || continue
      rhs -= L[k + 1] * result[order - k + 1]
    end
    result[order + 1] = apply_solve_plan(plan, rhs)
  end
  return result
end

function gram_feshbach_dressing(
  active_factor::MatrixSeries,
  cross::MatrixSeries,
  dark::MatrixSeries,
  N::Int,
  conditions::CompletionConditions,
)
  # If A = G G† and G Y = X, then X† A⁻¹ X = Y†Y. The same solve therefore
  # supplies both the active-to-dark dressing and the Feshbach/Schur residual.
  solved = lower_triangular_series_solve(active_factor, cross, N, conditions)
  correction = series_mul(series_adjoint(solved), solved, N)
  residual = series_sub(dark, correction, N)
  return residual, series_adjoint(solved)
end

function planned_undo_congruence_factor(
  factor::MatrixSeries, T::CompletionMatrix, N::Int, conditions::CompletionConditions
)
  validate_series_order(N)
  rows, columns = validate_matrix_series(factor)
  size(T, 1) == rows && size(T, 2) == rows ||
    throw(DimensionMismatch("congruence transform has incompatible dimensions"))

  adjoint_T = Matrix{CompletionScalar}(adjoint(T))
  plan = matrix_solve_plan(adjoint_T, conditions)
  result = Vector{CompletionMatrix}(undef, N + 1)
  for order in 0:N
    coefficient = matrix_coefficient(factor, order, rows, columns)
    result[order + 1] = apply_solve_plan(plan, coefficient)
  end
  return result
end

function swap_hermitian_coordinates!(
  A::CompletionMatrix, T::CompletionMatrix, first_index::Int, second_index::Int
)
  first_index == second_index && return nothing
  n = size(A, 1)
  for column in 1:n
    A[first_index, column], A[second_index, column] = A[second_index, column],
    A[first_index, column]
  end
  for row in 1:n
    A[row, first_index], A[row, second_index] = A[row, second_index], A[row, first_index]
  end
  for row in axes(T, 1)
    T[row, first_index], T[row, second_index] = T[row, second_index], T[row, first_index]
  end
  return nothing
end

function hermitian_congruence_eliminate_step!(
  reduced::CompletionMatrix, transform::CompletionMatrix, k::Int, pivot::CompletionScalar
)
  n = size(reduced, 1)
  coefficients = Vector{CompletionScalar}(undef, n - k)
  for column in (k + 1):n
    coefficients[column - k] = simplify_scalar(-reduced[k, column] / pivot)
  end

  # The trailing Hermitian block is the Schur update associated with the elementary
  # congruence. Update one triangle and restore the other by conjugation so no dense
  # elementary matrix multiplication is materialized.
  for column in (k + 1):n, row in column:n
    value = simplify_scalar(
      reduced[row, column] - reduced[row, k] * reduced[k, column] / pivot
    )
    if row == column
      reduced[row, column] = hermitian_real(value)
    else
      reduced[row, column] = value
      reduced[column, row] = simplify_scalar(conj(value))
    end
  end

  for column in (k + 1):n
    reduced[k, column] = completion_zero()
    reduced[column, k] = completion_zero()
  end

  # T <- T E, where E differs from the identity only in row k.
  for column in (k + 1):n
    coefficient = coefficients[column - k]
    for row in axes(transform, 1)
      transform[row, column] = simplify_scalar(
        transform[row, column] + coefficient * transform[row, k]
      )
    end
  end
  return nothing
end

function hermitian_eliminate_structured(
  A::CompletionMatrix, conditions::CompletionConditions
)
  hermitian_matrix(A) ||
    throw(ArgumentError("Hermitian elimination requires a Hermitian matrix"))
  n = size(A, 1)
  reduced = copy(A)
  transform = completion_identity(n)

  obstruction = zero_diagonal_obstruction(reduced, conditions)
  if obstruction != (0, 0)
    i, j = obstruction
    return HermitianElimination(
      transform,
      reduced,
      0,
      HERMITIAN_ZERO_DIAGONAL_COUPLING,
      simplify_scalar(reduced[i, j]),
    )
  end

  active_rank = 0
  for k in 1:n
    pivot_index = 0
    pivot_sign = SIGN_ZERO
    for candidate in k:n
      sign = structural_sign(reduced[candidate, candidate], conditions)
      if sign == SIGN_NEGATIVE
        return HermitianElimination(
          transform,
          reduced,
          active_rank,
          HERMITIAN_NEGATIVE_PIVOT,
          hermitian_real(reduced[candidate, candidate]),
        )
      elseif sign == SIGN_NONPOSITIVE
        return HermitianElimination(
          transform,
          reduced,
          active_rank,
          HERMITIAN_NONPOSITIVE_PIVOT,
          hermitian_real(reduced[candidate, candidate]),
        )
      elseif sign != SIGN_ZERO
        pivot_index = candidate
        pivot_sign = sign
        break
      end
    end

    pivot_index == 0 && break
    swap_hermitian_coordinates!(reduced, transform, k, pivot_index)
    pivot = hermitian_real(reduced[k, k])

    if pivot_sign == SIGN_UNKNOWN
      require_positivity!(conditions, pivot)
      require_regularity!(conditions, pivot)
    elseif pivot_sign == SIGN_NONNEGATIVE
      require_regularity!(conditions, pivot)
    end

    hermitian_congruence_eliminate_step!(reduced, transform, k, pivot)
    active_rank += 1

    obstruction = zero_diagonal_obstruction(reduced, conditions)
    if obstruction != (0, 0)
      i, j = obstruction
      return HermitianElimination(
        transform,
        reduced,
        active_rank,
        HERMITIAN_ZERO_DIAGONAL_COUPLING,
        simplify_scalar(reduced[i, j]),
      )
    end
  end

  return HermitianElimination(
    transform, reduced, active_rank, HERMITIAN_ELIMINATION_OK, completion_zero()
  )
end

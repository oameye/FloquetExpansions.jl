const CompletionScalar = Complex{Symbolics.Num}
const CompletionMatrix = Matrix{CompletionScalar}
const ScalarSeries = Vector{CompletionScalar}
const MatrixSeries = Vector{CompletionMatrix}

@enum StructuralSign::UInt8 begin
  SIGN_ZERO = 0x00
  SIGN_POSITIVE = 0x01
  SIGN_NEGATIVE = 0x02
  SIGN_NONNEGATIVE = 0x03
  SIGN_NONPOSITIVE = 0x04
  SIGN_UNKNOWN = 0x05
end

struct CompletionConditions
  positivity::Vector{CompletionScalar}
  regularity::Vector{CompletionScalar}
end

CompletionConditions() = CompletionConditions(CompletionScalar[], CompletionScalar[])

@enum HermitianEliminationStatus::UInt8 begin
  HERMITIAN_ELIMINATION_OK = 0x00
  HERMITIAN_NEGATIVE_PIVOT = 0x01
  HERMITIAN_NONPOSITIVE_PIVOT = 0x02
  HERMITIAN_ZERO_DIAGONAL_COUPLING = 0x03
end

struct HermitianElimination
  transform::CompletionMatrix
  reduced::CompletionMatrix
  active_rank::Int
  status::HermitianEliminationStatus
  obstruction::CompletionScalar
end

struct HermitianSeriesLDL
  lower::MatrixSeries
  diagonal::Vector{ScalarSeries}
end

completion_num(x::Symbolics.Num) = x
completion_num(x::Real) = Symbolics.Num(x)

completion_scalar(x::CompletionScalar) = x
completion_scalar(x::Symbolics.Num) = complex(x, completion_num(0))
completion_scalar(x::Real) = complex(completion_num(x), completion_num(0))
completion_scalar(x::Complex) = complex(completion_num(real(x)), completion_num(imag(x)))

completion_zero() = completion_scalar(0)
completion_one() = completion_scalar(1)

function completion_matrix_zeros(m::Int, n::Int)
  m >= 0 || throw(ArgumentError("matrix row count must be nonnegative"))
  n >= 0 || throw(ArgumentError("matrix column count must be nonnegative"))
  return fill(completion_zero(), m, n)
end

function completion_identity(n::Int)
  n >= 0 || throw(ArgumentError("matrix dimension must be nonnegative"))
  result = completion_matrix_zeros(n, n)
  for i in 1:n
    result[i, i] = completion_one()
  end
  return result
end

function simplify_scalar(z::CompletionScalar)
  return complex(
    completion_num(Symbolics.simplify(real(z))), completion_num(Symbolics.simplify(imag(z)))
  )
end

structurally_zero(x::Real) = iszero(x)
structurally_zero(x::Symbolics.Num) = issymzero(Symbolics.simplify(x))

function structurally_zero(z::CompletionScalar)
  simplified = simplify_scalar(z)
  return issymzero(real(simplified)) && issymzero(imag(simplified))
end

function structurally_zero(A::CompletionMatrix)
  return all(structurally_zero, A)
end

function structurally_equal(a::CompletionScalar, b::CompletionScalar)
  return structurally_zero(simplify_scalar(a - b))
end

function structurally_equal(A::CompletionMatrix, B::CompletionMatrix)
  size(A) == size(B) || return false
  for j in axes(A, 2), i in axes(A, 1)
    structurally_equal(A[i, j], B[i, j]) || return false
  end
  return true
end

function hermitian_real(z::CompletionScalar)
  simplified = simplify_scalar(z)
  issymzero(imag(simplified)) || throw(
    ArgumentError("Hermitian scalar must have structurally zero imaginary part; got `$z`")
  )
  return complex(real(simplified), completion_num(0))
end

function known_numeric_value(z::CompletionScalar)
  simplified = simplify_scalar(z)
  real_value = Symbolics.value(real(simplified))
  imag_value = Symbolics.value(imag(simplified))
  if real_value isa Real && imag_value isa Real
    return complex(real_value, imag_value)
  end
  return nothing
end

function condition_contains(conditions::Vector{CompletionScalar}, x::CompletionScalar)
  return any(p -> structurally_equal(p, x), conditions)
end

function require_positivity!(conditions::CompletionConditions, x::CompletionScalar)
  real_x = hermitian_real(x)
  condition_contains(conditions.positivity, real_x) || push!(conditions.positivity, real_x)
  return conditions
end

function require_regularity!(conditions::CompletionConditions, x::CompletionScalar)
  structurally_zero(x) &&
    throw(ArgumentError("cannot require a structurally zero pivot to be nonzero"))
  condition_contains(conditions.regularity, x) ||
    push!(conditions.regularity, simplify_scalar(x))
  return conditions
end

function structurally_nonzero(x::CompletionScalar, conditions::CompletionConditions)
  value = known_numeric_value(x)
  value !== nothing && return !iszero(value)
  return condition_contains(conditions.regularity, x)
end

function structural_sign(x::CompletionScalar, conditions::CompletionConditions)
  real_x = hermitian_real(x)
  structurally_zero(real_x) && return SIGN_ZERO

  value = known_numeric_value(real_x)
  if value !== nothing
    real_value = real(value)
    real_value > 0 && return SIGN_POSITIVE
    real_value < 0 && return SIGN_NEGATIVE
    return SIGN_ZERO
  end

  for p in conditions.positivity
    if structurally_equal(real_x, p)
      return condition_contains(conditions.regularity, p) ? SIGN_POSITIVE : SIGN_NONNEGATIVE
    elseif structurally_equal(real_x, -p)
      return condition_contains(conditions.regularity, p) ? SIGN_NEGATIVE : SIGN_NONPOSITIVE
    end
  end

  return SIGN_UNKNOWN
end

function validate_series_order(N::Int)
  N >= 0 || throw(ArgumentError("series order must be nonnegative"))
  return N
end

function validate_matrix_series(series::MatrixSeries)
  isempty(series) &&
    throw(ArgumentError("matrix series must contain at least one coefficient"))
  dims = size(first(series))
  all(size(A) == dims for A in series) ||
    throw(DimensionMismatch("all matrix-series coefficients must have the same dimensions"))
  return dims
end

function scalar_coefficient(series::ScalarSeries, n::Int)
  return n + 1 <= length(series) ? series[n + 1] : completion_zero()
end

function matrix_coefficient(series::MatrixSeries, n::Int, m::Int, p::Int)
  return n + 1 <= length(series) ? series[n + 1] : completion_matrix_zeros(m, p)
end

function series_add(a::ScalarSeries, b::ScalarSeries, N::Int)
  validate_series_order(N)
  result = Vector{CompletionScalar}(undef, N + 1)
  for n in 0:N
    result[n + 1] = simplify_scalar(scalar_coefficient(a, n) + scalar_coefficient(b, n))
  end
  return result
end

function series_sub(a::ScalarSeries, b::ScalarSeries, N::Int)
  validate_series_order(N)
  result = Vector{CompletionScalar}(undef, N + 1)
  for n in 0:N
    result[n + 1] = simplify_scalar(scalar_coefficient(a, n) - scalar_coefficient(b, n))
  end
  return result
end

function series_mul(a::ScalarSeries, b::ScalarSeries, N::Int)
  validate_series_order(N)
  result = fill(completion_zero(), N + 1)
  for n in 0:N
    coefficient = completion_zero()
    for k in 0:n
      coefficient += scalar_coefficient(a, k) * scalar_coefficient(b, n - k)
    end
    result[n + 1] = simplify_scalar(coefficient)
  end
  return result
end

function series_adjoint(a::ScalarSeries)
  return [simplify_scalar(conj(x)) for x in a]
end

function series_add(a::MatrixSeries, b::MatrixSeries, N::Int)
  validate_series_order(N)
  dims_a = validate_matrix_series(a)
  dims_b = validate_matrix_series(b)
  dims_a == dims_b || throw(DimensionMismatch("matrix-series dimensions must match"))
  m, p = dims_a
  result = Vector{CompletionMatrix}(undef, N + 1)
  for n in 0:N
    result[n + 1] = matrix_coefficient(a, n, m, p) + matrix_coefficient(b, n, m, p)
  end
  return result
end

function series_sub(a::MatrixSeries, b::MatrixSeries, N::Int)
  validate_series_order(N)
  dims_a = validate_matrix_series(a)
  dims_b = validate_matrix_series(b)
  dims_a == dims_b || throw(DimensionMismatch("matrix-series dimensions must match"))
  m, p = dims_a
  result = Vector{CompletionMatrix}(undef, N + 1)
  for n in 0:N
    result[n + 1] = matrix_coefficient(a, n, m, p) - matrix_coefficient(b, n, m, p)
  end
  return result
end

function series_mul(a::MatrixSeries, b::MatrixSeries, N::Int)
  validate_series_order(N)
  m, k_a = validate_matrix_series(a)
  k_b, p = validate_matrix_series(b)
  k_a == k_b || throw(DimensionMismatch("matrix-series inner dimensions must match"))
  result = [completion_matrix_zeros(m, p) for _ in 0:N]
  for n in 0:N
    coefficient = completion_matrix_zeros(m, p)
    for k in 0:n
      if k + 1 <= length(a) && n - k + 1 <= length(b)
        coefficient += a[k + 1] * b[n - k + 1]
      end
    end
    result[n + 1] = coefficient
  end
  return result
end

function series_adjoint(a::MatrixSeries)
  validate_matrix_series(a)
  return [Matrix{CompletionScalar}(adjoint(A)) for A in a]
end

function series_scale_shift(a::ScalarSeries, shift::Int, N::Int)
  validate_series_order(N)
  shift >= 0 || throw(ArgumentError("series shift must be nonnegative"))
  result = fill(completion_zero(), N + 1)
  for n in 0:min(N - shift, length(a) - 1)
    result[n + shift + 1] = a[n + 1]
  end
  return result
end

function series_scale_shift(a::MatrixSeries, shift::Int, N::Int)
  validate_series_order(N)
  shift >= 0 || throw(ArgumentError("series shift must be nonnegative"))
  m, p = validate_matrix_series(a)
  result = [completion_matrix_zeros(m, p) for _ in 0:N]
  for n in 0:min(N - shift, length(a) - 1)
    result[n + shift + 1] = copy(a[n + 1])
  end
  return result
end

function series_shift_down(a::ScalarSeries, shift::Int, N::Int)
  validate_series_order(N)
  shift >= 0 || throw(ArgumentError("series shift must be nonnegative"))
  for n in 0:min(shift - 1, length(a) - 1)
    structurally_zero(a[n + 1]) || throw(
      ArgumentError("cannot shift down a series with nonzero lower-order coefficients")
    )
  end
  result = fill(completion_zero(), N + 1)
  for n in 0:N
    n + shift + 1 <= length(a) || break
    result[n + 1] = a[n + shift + 1]
  end
  return result
end

function series_shift_down(a::MatrixSeries, shift::Int, N::Int)
  validate_series_order(N)
  shift >= 0 || throw(ArgumentError("series shift must be nonnegative"))
  m, p = validate_matrix_series(a)
  for n in 0:min(shift - 1, length(a) - 1)
    structurally_zero(a[n + 1]) || throw(
      ArgumentError("cannot shift down a series with nonzero lower-order coefficients")
    )
  end
  result = [completion_matrix_zeros(m, p) for _ in 0:N]
  for n in 0:N
    n + shift + 1 <= length(a) || break
    result[n + 1] = copy(a[n + shift + 1])
  end
  return result
end

function series_onset(a::ScalarSeries, N::Int)
  validate_series_order(N)
  for n in 0:min(N, length(a) - 1)
    structurally_zero(a[n + 1]) || return n
  end
  return -1
end

function scalar_series_inverse(a::ScalarSeries, N::Int, conditions::CompletionConditions)
  validate_series_order(N)
  isempty(a) && throw(ArgumentError("scalar series must contain at least one coefficient"))
  a0 = simplify_scalar(a[1])
  structurally_zero(a0) &&
    throw(ArgumentError("scalar series has a zero leading coefficient"))
  structurally_nonzero(a0, conditions) || require_regularity!(conditions, a0)

  result = fill(completion_zero(), N + 1)
  result[1] = simplify_scalar(completion_one() / a0)
  for n in 1:N
    convolution = completion_zero()
    for k in 1:n
      convolution += scalar_coefficient(a, k) * result[n - k + 1]
    end
    result[n + 1] = simplify_scalar(-result[1] * convolution)
  end
  return result
end

function scalar_series_sqrt(a::ScalarSeries, N::Int, conditions::CompletionConditions)
  validate_series_order(N)
  isempty(a) && throw(ArgumentError("scalar series must contain at least one coefficient"))
  real_series = [hermitian_real(x) for x in a]
  a0 = real_series[1]
  sign = structural_sign(a0, conditions)

  if sign == SIGN_NEGATIVE || sign == SIGN_NONPOSITIVE
    throw(
      ArgumentError("scalar square root requires a positive leading coefficient; got `$a0`")
    )
  elseif sign == SIGN_ZERO
    throw(ArgumentError("scalar square root requires a nonzero leading coefficient"))
  elseif sign == SIGN_NONNEGATIVE
    require_regularity!(conditions, a0)
  elseif sign == SIGN_UNKNOWN
    require_positivity!(conditions, a0)
    require_regularity!(conditions, a0)
  end

  result = fill(completion_zero(), N + 1)
  result[1] = completion_scalar(sqrt(real(a0)))
  for n in 1:N
    numerator = scalar_coefficient(real_series, n)
    for k in 1:(n - 1)
      numerator -= result[k + 1] * result[n - k + 1]
    end
    result[n + 1] = simplify_scalar(numerator / (2 * result[1]))
  end
  return result
end

function choose_pivot_row(
  A::CompletionMatrix, column::Int, conditions::CompletionConditions
)
  n = size(A, 1)
  fallback = 0
  for row in column:n
    pivot = A[row, column]
    structurally_zero(pivot) && continue
    structurally_nonzero(pivot, conditions) && return row
    fallback == 0 && (fallback = row)
  end
  return fallback
end

function matrix_solve(
  A::CompletionMatrix, B::CompletionMatrix, conditions::CompletionConditions
)
  n, m = size(A)
  n == m || throw(DimensionMismatch("solve matrix must be square"))
  size(B, 1) == n || throw(DimensionMismatch("right-hand side has incompatible row count"))

  rhs_columns = size(B, 2)
  work_A = copy(A)
  work_B = copy(B)

  for k in 1:n
    pivot_row = choose_pivot_row(work_A, k, conditions)
    pivot_row == 0 &&
      throw(ArgumentError("matrix is singular on the current symbolic stratum"))

    if pivot_row != k
      work_A[k, :], work_A[pivot_row, :] = copy(work_A[pivot_row, :]), copy(work_A[k, :])
      work_B[k, :], work_B[pivot_row, :] = copy(work_B[pivot_row, :]), copy(work_B[k, :])
    end

    pivot = simplify_scalar(work_A[k, k])
    structurally_nonzero(pivot, conditions) || require_regularity!(conditions, pivot)

    for i in (k + 1):n
      structurally_zero(work_A[i, k]) && continue
      factor = simplify_scalar(work_A[i, k] / pivot)
      work_A[i, k] = completion_zero()
      for j in (k + 1):n
        work_A[i, j] = simplify_scalar(work_A[i, j] - factor * work_A[k, j])
      end
      for j in 1:rhs_columns
        work_B[i, j] = simplify_scalar(work_B[i, j] - factor * work_B[k, j])
      end
    end
  end

  result = completion_matrix_zeros(n, rhs_columns)
  for column in 1:rhs_columns
    for i in n:-1:1
      residual = work_B[i, column]
      for j in (i + 1):n
        residual -= work_A[i, j] * result[j, column]
      end
      pivot = simplify_scalar(work_A[i, i])
      structurally_nonzero(pivot, conditions) || require_regularity!(conditions, pivot)
      result[i, column] = simplify_scalar(residual / pivot)
    end
  end
  return result
end

function lower_triangular_solve(
  L::CompletionMatrix, B::CompletionMatrix, conditions::CompletionConditions
)
  n, m = size(L)
  n == m || throw(DimensionMismatch("triangular solve matrix must be square"))
  size(B, 1) == n || throw(DimensionMismatch("right-hand side has incompatible row count"))
  result = completion_matrix_zeros(n, size(B, 2))
  for column in axes(B, 2), i in 1:n
    residual = B[i, column]
    for j in 1:(i - 1)
      residual -= L[i, j] * result[j, column]
    end
    pivot = simplify_scalar(L[i, i])
    structurally_nonzero(pivot, conditions) || require_regularity!(conditions, pivot)
    result[i, column] = simplify_scalar(residual / pivot)
  end
  return result
end

function unit_lower_triangular_solve(L::CompletionMatrix, B::CompletionMatrix)
  n, m = size(L)
  n == m || throw(DimensionMismatch("triangular solve matrix must be square"))
  size(B, 1) == n || throw(DimensionMismatch("right-hand side has incompatible row count"))
  result = completion_matrix_zeros(n, size(B, 2))
  for column in axes(B, 2), i in 1:n
    residual = B[i, column]
    for j in 1:(i - 1)
      residual -= L[i, j] * result[j, column]
    end
    result[i, column] = simplify_scalar(residual)
  end
  return result
end

function adjoint_triangular_solve(
  L::CompletionMatrix, B::CompletionMatrix, conditions::CompletionConditions
)
  n, m = size(L)
  n == m || throw(DimensionMismatch("triangular solve matrix must be square"))
  size(B, 1) == n || throw(DimensionMismatch("right-hand side has incompatible row count"))
  result = completion_matrix_zeros(n, size(B, 2))
  for column in axes(B, 2), i in n:-1:1
    residual = B[i, column]
    for j in (i + 1):n
      residual -= conj(L[j, i]) * result[j, column]
    end
    pivot = simplify_scalar(conj(L[i, i]))
    structurally_nonzero(pivot, conditions) || require_regularity!(conditions, pivot)
    result[i, column] = simplify_scalar(residual / pivot)
  end
  return result
end

function series_solve(
  A::MatrixSeries, B::MatrixSeries, N::Int, conditions::CompletionConditions
)
  validate_series_order(N)
  n, m = validate_matrix_series(A)
  n == m || throw(DimensionMismatch("left matrix series must be square"))
  b_rows, b_columns = validate_matrix_series(B)
  b_rows == n ||
    throw(DimensionMismatch("matrix-series right-hand side has incompatible rows"))

  result = [completion_matrix_zeros(n, b_columns) for _ in 0:N]
  A0 = A[1]
  for order in 0:N
    rhs = copy(matrix_coefficient(B, order, b_rows, b_columns))
    for k in 1:order
      k + 1 <= length(A) || continue
      rhs -= A[k + 1] * result[order - k + 1]
    end
    result[order + 1] = matrix_solve(A0, rhs, conditions)
  end
  return result
end

function series_inverse(A::MatrixSeries, N::Int, conditions::CompletionConditions)
  n, m = validate_matrix_series(A)
  n == m || throw(DimensionMismatch("matrix series must be square"))
  identity_series = [completion_matrix_zeros(n, n) for _ in 0:N]
  identity_series[1] = completion_identity(n)
  return series_solve(A, identity_series, N, conditions)
end

function series_schur(
  A::MatrixSeries,
  X::MatrixSeries,
  C::MatrixSeries,
  N::Int,
  conditions::CompletionConditions,
)
  solved = series_solve(A, X, N, conditions)
  correction = series_mul(series_adjoint(X), solved, N)
  return series_sub(C, correction, N)
end

function apply_congruence(series::MatrixSeries, T::CompletionMatrix, N::Int)
  validate_series_order(N)
  n, m = validate_matrix_series(series)
  n == m ||
    throw(DimensionMismatch("congruence requires square matrix-series coefficients"))
  size(T) == (n, n) ||
    throw(DimensionMismatch("congruence transform must be square and match the series"))
  result = Vector{CompletionMatrix}(undef, N + 1)
  for order in 0:N
    coefficient = matrix_coefficient(series, order, n, n)
    result[order + 1] = Matrix{CompletionScalar}(adjoint(T) * coefficient * T)
  end
  return result
end

function undo_congruence_factor(
  factor::MatrixSeries, T::CompletionMatrix, N::Int, conditions::CompletionConditions
)
  validate_series_order(N)
  rows, columns = validate_matrix_series(factor)
  size(T, 1) == rows && size(T, 2) == rows ||
    throw(DimensionMismatch("congruence transform has incompatible dimensions"))
  adjoint_T = Matrix{CompletionScalar}(adjoint(T))
  result = Vector{CompletionMatrix}(undef, N + 1)
  for order in 0:N
    coefficient = matrix_coefficient(factor, order, rows, columns)
    result[order + 1] = matrix_solve(adjoint_T, coefficient, conditions)
  end
  return result
end

function hermitian_matrix(A::CompletionMatrix)
  size(A, 1) == size(A, 2) || return false
  for j in axes(A, 2), i in 1:j
    structurally_equal(A[i, j], conj(A[j, i])) || return false
  end
  return true
end

function hermitian_series(A::MatrixSeries)
  validate_matrix_series(A)
  return all(hermitian_matrix, A)
end

function zero_diagonal_obstruction(A::CompletionMatrix, conditions::CompletionConditions)
  hermitian_matrix(A) ||
    throw(ArgumentError("Hermitian elimination requires a Hermitian matrix"))
  n = size(A, 1)
  for i in 1:n
    structurally_zero(hermitian_real(A[i, i])) || continue
    for j in 1:n
      i == j && continue
      coupling = A[i, j]
      structurally_zero(coupling) && continue
      structurally_nonzero(coupling, conditions) ||
        require_regularity!(conditions, coupling)
      return (i, j)
    end
  end
  return (0, 0)
end

function swap_congruence_coordinates!(
  A::CompletionMatrix, T::CompletionMatrix, i::Int, j::Int
)
  i == j && return nothing
  A[i, :], A[j, :] = copy(A[j, :]), copy(A[i, :])
  A[:, i], A[:, j] = copy(A[:, j]), copy(A[:, i])
  T[:, i], T[:, j] = copy(T[:, j]), copy(T[:, i])
  return nothing
end

function hermitian_eliminate(A::CompletionMatrix, conditions::CompletionConditions)
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
    swap_congruence_coordinates!(reduced, transform, k, pivot_index)
    pivot = hermitian_real(reduced[k, k])

    if pivot_sign == SIGN_UNKNOWN
      require_positivity!(conditions, pivot)
      require_regularity!(conditions, pivot)
    elseif pivot_sign == SIGN_NONNEGATIVE
      require_regularity!(conditions, pivot)
    end

    E = completion_identity(n)
    for j in (k + 1):n
      structurally_zero(reduced[k, j]) && continue
      E[k, j] = simplify_scalar(-reduced[k, j] / pivot)
    end
    reduced = Matrix{CompletionScalar}(adjoint(E) * reduced * E)
    transform = transform * E
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

function scalar_series_entry(A::MatrixSeries, i::Int, j::Int, N::Int)
  rows, columns = validate_matrix_series(A)
  1 <= i <= rows || throw(BoundsError(first(A), (i, j)))
  1 <= j <= columns || throw(BoundsError(first(A), (i, j)))
  result = fill(completion_zero(), N + 1)
  for order in 0:N
    order + 1 <= length(A) || break
    result[order + 1] = A[order + 1][i, j]
  end
  return result
end

function set_scalar_series_entry!(
  A::MatrixSeries, i::Int, j::Int, values::ScalarSeries, N::Int
)
  for order in 0:N
    A[order + 1][i, j] = scalar_coefficient(values, order)
  end
  return A
end

function scaled_scalar_product(a::ScalarSeries, b::ScalarSeries, c::ScalarSeries, N::Int)
  return series_mul(series_mul(a, b, N), c, N)
end

function graded_ldl(A::MatrixSeries, N::Int, conditions::CompletionConditions)
  validate_series_order(N)
  rows, columns = validate_matrix_series(A)
  rows == columns ||
    throw(DimensionMismatch("graded LDL factorization requires square matrices"))
  hermitian_series(A) ||
    throw(ArgumentError("graded LDL factorization requires a Hermitian series"))

  lower = [completion_matrix_zeros(rows, rows) for _ in 0:N]
  for i in 1:rows
    lower[1][i, i] = completion_one()
  end
  diagonal = Vector{ScalarSeries}(undef, rows)

  for j in 1:rows
    delta = scalar_series_entry(A, j, j, N)
    for k in 1:(j - 1)
      Ljk = scalar_series_entry(lower, j, k, N)
      delta = series_sub(
        delta, scaled_scalar_product(Ljk, diagonal[k], series_adjoint(Ljk), N), N
      )
    end
    delta = [hermitian_real(x) for x in delta]

    sign = structural_sign(delta[1], conditions)
    if sign == SIGN_NEGATIVE || sign == SIGN_NONPOSITIVE
      throw(
        ArgumentError("graded LDL factorization encountered a negative Hermitian pivot")
      )
    elseif sign == SIGN_ZERO
      throw(ArgumentError("graded LDL factorization encountered a dark leading pivot"))
    elseif sign == SIGN_UNKNOWN
      require_positivity!(conditions, delta[1])
      require_regularity!(conditions, delta[1])
    elseif sign == SIGN_NONNEGATIVE
      require_regularity!(conditions, delta[1])
    end
    diagonal[j] = delta

    inverse_delta = scalar_series_inverse(delta, N, conditions)
    for i in (j + 1):rows
      numerator = scalar_series_entry(A, i, j, N)
      for k in 1:(j - 1)
        lower_ik = scalar_series_entry(lower, i, k, N)
        Ljk = scalar_series_entry(lower, j, k, N)
        numerator = series_sub(
          numerator, scaled_scalar_product(lower_ik, diagonal[k], series_adjoint(Ljk), N), N
        )
      end
      set_scalar_series_entry!(lower, i, j, series_mul(numerator, inverse_delta, N), N)
    end
  end

  return HermitianSeriesLDL(lower, diagonal)
end

function diagonal_matrix_series(diagonal::Vector{ScalarSeries}, N::Int)
  n = length(diagonal)
  result = [completion_matrix_zeros(n, n) for _ in 0:N]
  for j in 1:n, order in 0:N
    result[order + 1][j, j] = scalar_coefficient(diagonal[j], order)
  end
  return result
end

function reconstruct_ldl(factorization::HermitianSeriesLDL, N::Int)
  diagonal = diagonal_matrix_series(factorization.diagonal, N)
  return series_mul(
    series_mul(factorization.lower, diagonal, N), series_adjoint(factorization.lower), N
  )
end

function ldl_gram_factor(
  factorization::HermitianSeriesLDL, N::Int, conditions::CompletionConditions
)
  n = length(factorization.diagonal)
  square_root = [completion_matrix_zeros(n, n) for _ in 0:N]
  for j in 1:n
    root = scalar_series_sqrt(factorization.diagonal[j], N, conditions)
    for order in 0:N
      square_root[order + 1][j, j] = root[order + 1]
    end
  end
  return series_mul(factorization.lower, square_root, N)
end

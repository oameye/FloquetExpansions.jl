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

function completion_symbolic_zero(x::Symbolics.Num)::Bool
  simplified = Symbolics.simplify(x)::Symbolics.Num
  return isequal(simplified, completion_num(0))
end

structurally_zero(x::Symbolics.Num) = completion_symbolic_zero(x)

function structurally_zero(z::CompletionScalar)
  simplified = simplify_scalar(z)
  return completion_symbolic_zero(real(simplified)) &&
         completion_symbolic_zero(imag(simplified))
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
  completion_symbolic_zero(imag(simplified)) || throw(
    ArgumentError("Hermitian scalar must have structurally zero imaginary part; got `$z`")
  )
  return complex(real(simplified), completion_num(0))
end

symbolic_truth(value::Bool) = value

function symbolic_truth(value::Symbolics.Num)::Bool
  simplified = Symbolics.simplify(value)::Symbolics.Num
  return Symbolics.value(simplified) === true
end

function symbolically_positive(value::Symbolics.Num)::Bool
  return symbolic_truth(value > 0)
end

function symbolically_negative(value::Symbolics.Num)::Bool
  return symbolic_truth(value < 0)
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
  simplified = simplify_scalar(x)
  structurally_zero(simplified) && return false
  real_part = real(simplified)
  imag_part = imag(simplified)
  if symbolically_positive(real_part) ||
    symbolically_negative(real_part) ||
    symbolically_positive(imag_part) ||
    symbolically_negative(imag_part)
    return true
  end
  return condition_contains(conditions.regularity, simplified)
end

function structural_sign(x::CompletionScalar, conditions::CompletionConditions)
  real_x = hermitian_real(x)
  structurally_zero(real_x) && return SIGN_ZERO

  real_value = real(real_x)
  symbolically_positive(real_value) && return SIGN_POSITIVE
  symbolically_negative(real_value) && return SIGN_NEGATIVE

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

function graded_ldl_diagonal(
  A::MatrixSeries, lower::MatrixSeries, diagonal::Vector{ScalarSeries}, j::Int, N::Int
)::ScalarSeries
  delta = scalar_series_entry(A, j, j, N)
  for k in 1:(j - 1)
    lower_jk = scalar_series_entry(lower, j, k, N)
    delta = series_sub(
      delta, scaled_scalar_product(lower_jk, diagonal[k], series_adjoint(lower_jk), N), N
    )
  end
  for index in eachindex(delta)
    delta[index] = hermitian_real(delta[index])
  end
  return delta
end

function require_positive_ldl_pivot!(delta::ScalarSeries, conditions::CompletionConditions)
  sign = structural_sign(delta[1], conditions)
  if sign == SIGN_NEGATIVE || sign == SIGN_NONPOSITIVE
    throw(ArgumentError("graded LDL factorization encountered a negative Hermitian pivot"))
  elseif sign == SIGN_ZERO
    throw(ArgumentError("graded LDL factorization encountered a dark leading pivot"))
  elseif sign == SIGN_UNKNOWN
    require_positivity!(conditions, delta[1])
    require_regularity!(conditions, delta[1])
  elseif sign == SIGN_NONNEGATIVE
    require_regularity!(conditions, delta[1])
  end
  return conditions
end

function graded_ldl_numerator(
  A::MatrixSeries,
  lower::MatrixSeries,
  diagonal::Vector{ScalarSeries},
  i::Int,
  j::Int,
  N::Int,
)::ScalarSeries
  numerator = scalar_series_entry(A, i, j, N)
  for k in 1:(j - 1)
    lower_ik = scalar_series_entry(lower, i, k, N)
    lower_jk = scalar_series_entry(lower, j, k, N)
    numerator = series_sub(
      numerator,
      scaled_scalar_product(lower_ik, diagonal[k], series_adjoint(lower_jk), N),
      N,
    )
  end
  return numerator
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
    delta = graded_ldl_diagonal(A, lower, diagonal, j, N)
    require_positive_ldl_pivot!(delta, conditions)
    diagonal[j] = delta

    inverse_delta = scalar_series_inverse(delta, N, conditions)
    for i in (j + 1):rows
      numerator = graded_ldl_numerator(A, lower, diagonal, i, j, N)
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

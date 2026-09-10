from pathlib import Path

p = Path('src/matrix_series.jl')
s = p.read_text()
old = '''function series_mul(a::MatrixSeries, b::MatrixSeries, N::Int)
  validate_series_order(N)
  m, k_a = validate_matrix_series(a)
  k_b, p = validate_matrix_series(b)
  k_a == k_b || throw(DimensionMismatch("matrix-series inner dimensions must match"))
  result = Vector{CompletionMatrix}(undef, N + 1)
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
'''
new = '''function matrix_product_dimensions(
  target::CompletionMatrix, left::CompletionMatrix, right::CompletionMatrix
)
  rows, inner = size(left)
  right_inner, columns = size(right)
  inner == right_inner || throw(DimensionMismatch("matrix product inner dimensions must match"))
  size(target) == (rows, columns) ||
    throw(DimensionMismatch("matrix product target has incompatible dimensions"))
  return rows, inner, columns
end

function add_matrix_product!(
  target::CompletionMatrix, left::CompletionMatrix, right::CompletionMatrix
)
  rows, inner, columns = matrix_product_dimensions(target, left, right)
  for column in 1:columns, row in 1:rows
    product = completion_zero()
    for index in 1:inner
      product += left[row, index] * right[index, column]
    end
    target[row, column] += product
  end
  return target
end

function subtract_matrix_product!(
  target::CompletionMatrix, left::CompletionMatrix, right::CompletionMatrix
)
  rows, inner, columns = matrix_product_dimensions(target, left, right)
  for column in 1:columns, row in 1:rows
    product = completion_zero()
    for index in 1:inner
      product += left[row, index] * right[index, column]
    end
    target[row, column] -= product
  end
  return target
end

function series_mul(a::MatrixSeries, b::MatrixSeries, N::Int)
  validate_series_order(N)
  m, k_a = validate_matrix_series(a)
  k_b, p = validate_matrix_series(b)
  k_a == k_b || throw(DimensionMismatch("matrix-series inner dimensions must match"))
  result = Vector{CompletionMatrix}(undef, N + 1)
  for n in 0:N
    coefficient = completion_matrix_zeros(m, p)
    for k in 0:n
      if k + 1 <= length(a) && n - k + 1 <= length(b)
        add_matrix_product!(coefficient, a[k + 1], b[n - k + 1])
      end
    end
    result[n + 1] = coefficient
  end
  return result
end
'''
assert s.count(old) == 1, 'matrix series target not found exactly once'
p.write_text(s.replace(old, new))

p = Path('src/completion_linear_algebra.jl')
s = p.read_text()
old1 = '      rhs -= A[k + 1] * result[order - k + 1]\n'
old2 = '      rhs -= L[k + 1] * result[order - k + 1]\n'
assert s.count(old1) == 1, 'planned solve target not found exactly once'
assert s.count(old2) == 1, 'lower triangular solve target not found exactly once'
s = s.replace(old1, '      subtract_matrix_product!(rhs, A[k + 1], result[order - k + 1])\n')
s = s.replace(old2, '      subtract_matrix_product!(rhs, L[k + 1], result[order - k + 1])\n')
p.write_text(s)

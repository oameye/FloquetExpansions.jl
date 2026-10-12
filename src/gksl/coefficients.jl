const KossakowskiMatrix = Matrix{SQA.CNum}

@inline coefficient_zero() = convert(SQA.CNum, 0)
@inline coefficient_one() = convert(SQA.CNum, 1)
@inline simplify_coefficient(value::SQA.CNum) = SQA.simplify(value)

coefficient_matrix(rows::Int, columns::Int) = fill(coefficient_zero(), rows, columns)

function simplify_matrix!(matrix::KossakowskiMatrix)
  for index in eachindex(matrix)
    matrix[index] = simplify_coefficient(matrix[index])
  end
  return matrix
end

function multiply_coefficients(left::KossakowskiMatrix, right::KossakowskiMatrix)
  size(left, 2) == size(right, 1) ||
    throw(DimensionMismatch("matrix dimensions do not match"))
  result = coefficient_matrix(size(left, 1), size(right, 2))
  for row in axes(result, 1), column in axes(result, 2)
    value = coefficient_zero()
    for index in axes(left, 2)
      value = value + left[row, index] * right[index, column]
    end
    result[row, column] = simplify_coefficient(value)
  end
  return result
end

function adjoint_coefficients(matrix::KossakowskiMatrix)
  result = coefficient_matrix(size(matrix, 2), size(matrix, 1))
  for row in axes(matrix, 1), column in axes(matrix, 2)
    result[column, row] = conj(matrix[row, column])
  end
  return result
end

function conjugate_pair(left::SQA.CNum, right::SQA.CNum)
  difference = SQA.to_num(left) - conj(SQA.to_num(right))
  return iszero(Symbolics.expand(real(difference))) &&
         iszero(Symbolics.expand(imag(difference)))
end

function matrix_is_hermitian(matrix::KossakowskiMatrix)
  size(matrix, 1) == size(matrix, 2) || return false
  for row in axes(matrix, 1), column in row:size(matrix, 2)
    conjugate_pair(matrix[row, column], matrix[column, row]) && continue
    iszero(simplify_coefficient(matrix[row, column] - conj(matrix[column, row]))) ||
      return false
  end
  return true
end

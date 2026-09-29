# Internal Gram geometry shared by positive completion and generator-native GKLS work.
#
# These helpers operate on a fixed amplitude frame. They do not choose a Floquet/static gauge and
# they do not perform positive completion of an already-computed effective generator.

function gram_active_frame(
  B::AbstractMatrix{T}; rtol::Real=1e-10, atol::Real=0.0
) where {T<:Number}
  rows, columns = size(B)
  decomposition = LinearAlgebra.svd(Matrix(B); full=true)
  scale = isempty(decomposition.S) ? 0.0 : maximum(decomposition.S)
  tolerance = max(float(atol), float(rtol) * max(1.0, float(scale)))
  rank = count(>(tolerance), decomposition.S)

  U = decomposition.U
  V = Matrix(adjoint(decomposition.Vt))
  active = U[:, 1:rank]
  dark = U[:, (rank + 1):rows]
  channel = V[:, 1:rank]
  singular = decomposition.S[1:rank]
  return (; active, dark, channel, singular, rank, tolerance, rows, columns)
end

function gram_tangent_lift(
  B::AbstractMatrix{T}, target::AbstractMatrix{T}; rtol::Real=1e-10, atol::Real=0.0
) where {T<:Number}
  rows, columns = size(B)
  size(target) == (rows, rows) ||
    throw(DimensionMismatch("tangent Gram target must match the amplitude row dimension"))

  hermitian_target = (target + adjoint(target)) / 2
  scale = max(1.0, float(LinearAlgebra.norm(hermitian_target)))
  hermitian_error = LinearAlgebra.norm(target - adjoint(target))
  hermitian_error <= max(float(atol), float(rtol) * scale) ||
    throw(ArgumentError("tangent Gram target must be Hermitian"))

  frame = gram_active_frame(B; rtol, atol)
  if frame.rank == 0
    LinearAlgebra.norm(hermitian_target) <= max(float(atol), float(rtol) * scale) ||
      throw(ArgumentError("nonzero tangent Gram target with no active amplitude columns"))
    return zeros(T, rows, columns)
  end

  U = hcat(frame.active, frame.dark)
  transformed = adjoint(U) * hermitian_target * U
  rank = frame.rank
  sigma = frame.singular

  bright = transformed[1:rank, 1:rank]
  dark_bright = transformed[(rank + 1):rows, 1:rank]

  A = Matrix{T}(undef, rank, rank)
  for j in 1:rank, i in 1:rank
    denominator = sigma[i]^2 + sigma[j]^2
    A[i, j] = bright[i, j] * (sigma[i] * sigma[j] / denominator)
  end
  X = LinearAlgebra.Diagonal(inv.(sigma)) * A
  Z = if isempty(frame.dark)
    zeros(T, 0, rank)
  else
    dark_bright * LinearAlgebra.Diagonal(inv.(sigma))
  end
  reduced = vcat(X, Z)
  correction = U * reduced * adjoint(frame.channel)

  reconstructed = B * adjoint(correction) + correction * adjoint(B)
  residual = hermitian_target - reconstructed
  tolerance = max(float(atol), float(rtol) * scale)
  LinearAlgebra.norm(residual) <= 20 * tolerance ||
    throw(ArgumentError("target is outside the tangent space of the active Gram factor"))

  return correction
end

function positive_gram_factor(
  form::AbstractMatrix{T}; rtol::Real=1e-10, atol::Real=0.0
) where {T<:Number}
  rows, columns = size(form)
  rows == columns || throw(DimensionMismatch("positive Gram form must be square"))

  hermitian_form = (form + adjoint(form)) / 2
  scale = max(1.0, float(LinearAlgebra.norm(hermitian_form)))
  tolerance = max(float(atol), float(rtol) * scale)
  LinearAlgebra.norm(form - adjoint(form)) <= tolerance ||
    throw(ArgumentError("positive Gram form must be Hermitian"))

  decomposition = LinearAlgebra.eigen(LinearAlgebra.Hermitian(hermitian_form))
  minimum_value = isempty(decomposition.values) ? 0.0 : minimum(decomposition.values)
  minimum_value >= -tolerance ||
    throw(ArgumentError("Gram form is not positive semidefinite on the resolved stratum"))

  keep = findall(>(tolerance), decomposition.values)
  isempty(keep) && return zeros(T, rows, 0)
  values = max.(decomposition.values[keep], zero(eltype(decomposition.values)))
  return decomposition.vectors[:, keep] * LinearAlgebra.Diagonal(sqrt.(values))
end

function gram_feshbach_residual(dark::MatrixSeries, solved::MatrixSeries, N::Int)
  correction = series_mul(series_adjoint(solved), solved, N)
  return series_sub(dark, correction, N)
end

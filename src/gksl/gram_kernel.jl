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

function gram_hermitian_coordinates(form::AbstractMatrix)
  rows, columns = size(form)
  rows == columns || throw(DimensionMismatch("Hermitian coordinate input must be square"))
  values = Float64[]
  for i in 1:rows
    push!(values, real(form[i, i]))
  end
  for i in 1:rows, j in (i + 1):rows
    push!(values, sqrt(2) * real(form[i, j]))
    push!(values, sqrt(2) * imag(form[i, j]))
  end
  return values
end

function gram_hermitian_from_coordinates(values::AbstractVector{<:Real}, n::Int)
  length(values) == n^2 ||
    throw(DimensionMismatch("Hermitian coordinate vector must have length n^2"))
  form = zeros(ComplexF64, n, n)
  k = 1
  for i in 1:n
    form[i, i] = values[k]
    k += 1
  end
  for i in 1:n, j in (i + 1):n
    form[i, j] = (values[k] + im * values[k + 1]) / sqrt(2)
    form[j, i] = conj(form[i, j])
    k += 2
  end
  return form
end

function gram_psd_projection(values::AbstractVector{<:Real}, n::Int)
  form = gram_hermitian_from_coordinates(values, n)
  decomposition = LinearAlgebra.eigen(LinearAlgebra.Hermitian(form))
  projected_values = max.(decomposition.values, zero(eltype(decomposition.values)))
  projected =
    decomposition.vectors *
    LinearAlgebra.Diagonal(projected_values) *
    adjoint(decomposition.vectors)
  return gram_hermitian_coordinates(projected)
end

function gram_is_psd(form::AbstractMatrix, tolerance::Real)
  values = LinearAlgebra.eigvals(LinearAlgebra.Hermitian((form + adjoint(form)) / 2))
  isempty(values) && return true
  scale = max(1.0, float(LinearAlgebra.norm(form)))
  return minimum(values) >= -float(tolerance) * scale
end

function gram_affine_geometry(
  delta::AbstractMatrix, Phi::AbstractMatrix{<:Real}, rtol::Real
)
  rows, columns = size(delta)
  rows == columns || throw(DimensionMismatch("affine PSD residual must be square"))

  hermitian_delta = (delta + adjoint(delta)) / 2
  scale = max(1.0, float(LinearAlgebra.norm(hermitian_delta)))
  LinearAlgebra.norm(delta - adjoint(delta)) <= float(rtol) * scale ||
    throw(ArgumentError("affine PSD residual must be Hermitian"))

  delta_coordinates = gram_hermitian_coordinates(hermitian_delta)
  size(Phi, 1) == length(delta_coordinates) || throw(
    DimensionMismatch("affine PSD map row count must match Hermitian coordinate dimension"),
  )

  Phi_matrix = Matrix{Float64}(Phi)
  Phi_pinv = if size(Phi_matrix, 2) == 0
    zeros(Float64, 0, length(delta_coordinates))
  else
    LinearAlgebra.pinv(Phi_matrix; rtol=float(rtol))
  end
  return (; rows, delta_coordinates, Phi_matrix, Phi_pinv)
end

function gram_affine_projection(values, geometry)
  delta = geometry.delta_coordinates
  return delta + geometry.Phi_matrix * (geometry.Phi_pinv * (values - delta))
end

function gram_affine_result(vector, geometry, canonical, iterations)
  coordinates = geometry.Phi_pinv * (vector - geometry.delta_coordinates)
  reconstructed = geometry.delta_coordinates + geometry.Phi_matrix * coordinates
  form = gram_hermitian_from_coordinates(reconstructed, geometry.rows)
  return (; coordinates, form, canonical, iterations)
end

function gram_canonical_affine_result(geometry, tolerance)
  coordinates = -geometry.Phi_pinv * geometry.delta_coordinates
  vector = geometry.delta_coordinates + geometry.Phi_matrix * coordinates
  form = gram_hermitian_from_coordinates(vector, geometry.rows)
  gram_is_psd(form, tolerance) || return nothing
  return (; coordinates, form, canonical=true, iterations=0)
end

function gram_dykstra_affine_psd(geometry, tolerance, maxiter)
  current = zeros(Float64, length(geometry.delta_coordinates))
  affine_correction = zeros(Float64, length(current))
  positive_correction = zeros(Float64, length(current))

  for iteration in 1:maxiter
    shifted_affine = current + affine_correction
    affine = gram_affine_projection(shifted_affine, geometry)
    affine_correction = shifted_affine - affine

    shifted_positive = affine + positive_correction
    positive = gram_psd_projection(shifted_positive, geometry.rows)
    positive_correction = shifted_positive - positive
    current = positive

    candidate = gram_affine_projection(current, geometry)
    scale = max(1.0, LinearAlgebra.norm(candidate), LinearAlgebra.norm(current))
    distance = LinearAlgebra.norm(candidate - current)
    form = gram_hermitian_from_coordinates(candidate, geometry.rows)
    distance <= 20 * tolerance * scale || continue
    gram_is_psd(form, tolerance) || continue
    return gram_affine_result(candidate, geometry, false, iteration)
  end

  return throw(
    ArgumentError("affine Hermitian slice did not reach the PSD cone within tolerance")
  )
end

function positive_affine_section(
  delta::AbstractMatrix, Phi::AbstractMatrix{<:Real}; rtol::Real=1e-10, maxiter::Int=10_000
)
  maxiter > 0 || throw(ArgumentError("affine PSD solver requires maxiter > 0"))
  geometry = gram_affine_geometry(delta, Phi, rtol)
  tolerance = float(rtol)
  canonical = gram_canonical_affine_result(geometry, tolerance)
  canonical === nothing || return canonical
  return gram_dykstra_affine_psd(geometry, tolerance, maxiter)
end

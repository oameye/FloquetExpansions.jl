struct ExactStaticSolution{T,R}
  coordinates::Vector{R}
  newborn::Matrix{T}
  newborn_weights::Vector{R}
  correction::Matrix{T}
  coefficient::Matrix{T}
  dark_residual::Matrix{T}
end

exact_zeros(::Type{S}, dims::Int...) where {S} = fill!(Array{S}(undef, dims...), zero(S))

function exact_solve(A::AbstractMatrix{S}, B::AbstractVecOrMat{S}) where {S}
  n = size(A, 1)
  size(A, 2) == n || throw(DimensionMismatch("exact solve needs a square matrix"))
  M = Matrix{S}(A)
  X = Array{S}(B)
  for column in 1:n
    offset = findfirst(!iszero, view(M, column:n, column))
    offset === nothing && throw(ArgumentError("exact solve met a singular matrix"))
    pivot = column + offset - 1
    if pivot != column
      M[[column, pivot], :] = M[[pivot, column], :]
      X[[column, pivot], :] = X[[pivot, column], :]
    end
    scale = inv(M[column, column])
    M[column, :] *= scale
    X[column, :] *= scale
    for r in 1:n
      r == column && continue
      factor = M[r, column]
      iszero(factor) && continue
      M[r, :] -= factor * M[column, :]
      X[r, :] -= factor * X[column, :]
    end
  end
  return X
end

function exact_row_basis(A::AbstractMatrix{T}) where {T}
  basis = Vector{T}[]
  pivots = Int[]
  selected = Int[]
  for r in axes(A, 1)
    v = Vector{T}(A[r, :])
    for (b, p) in zip(basis, pivots)
      iszero(v[p]) || (v -= (v[p] / b[p]) * b)
    end
    p = findfirst(!iszero, v)
    p === nothing && continue
    push!(basis, v)
    push!(pivots, p)
    push!(selected, r)
  end
  return selected
end

function exact_min_norm_solve(
  A::AbstractMatrix{R}, b::AbstractVector{R}, metric::AbstractMatrix{R}
) where {R}
  rows = exact_row_basis(A)
  isempty(rows) && return exact_zeros(R, size(A, 2))
  Ar = A[rows, :]
  weighted = exact_solve(metric, Matrix(transpose(Ar)))
  solution = weighted * exact_solve(Ar * weighted, b[rows])
  A * solution == b || throw(ArgumentError("exact linear system is inconsistent"))
  return solution
end

function exact_hermitian_units(::Type{T}, n::Int) where {T}
  units = Matrix{T}[]
  for i in 1:n
    E = exact_zeros(T, n, n)
    E[i, i] = 1
    push!(units, E)
  end
  for i in 1:n, j in (i + 1):n
    E = exact_zeros(T, n, n)
    E[i, j] = 1
    E[j, i] = 1
    push!(units, E)
    F = exact_zeros(T, n, n)
    F[i, j] = im
    F[j, i] = -im
    push!(units, F)
  end
  return units
end

function exact_hermitian_coordinates(X::AbstractMatrix{T}, ::Type{R}) where {T,R}
  n = size(X, 1)
  values = R[real(X[i, i]) for i in 1:n]
  for i in 1:n, j in (i + 1):n
    push!(values, real(X[i, j]), imag(X[i, j]))
  end
  return values
end

function exact_form_metric(
  G::AbstractMatrix{T}, units::Vector{Matrix{T}}, ::Type{R}
) where {T,R}
  n = length(units)
  supports = [Tuple.(findall(!iszero, E)) for E in units]
  M = exact_zeros(R, n, n)
  for k in 1:n, l in k:n
    value = zero(T)
    for (b, c) in supports[k], (e, a) in supports[l]
      value += G[a, b] * units[k][b, c] * G[c, e] * units[l][e, a]
    end
    M[k, l] = real(value)
    M[l, k] = real(value)
  end
  return M
end

function exact_dark_projector(B::AbstractMatrix{T}, G::AbstractMatrix{T}) where {T}
  n = size(G, 1)
  size(B, 2) == 0 && return Matrix{T}(LinearAlgebra.I, n, n)
  return Matrix{T}(LinearAlgebra.I, n, n) -
         B * exact_solve(adjoint(B) * G * B, adjoint(B) * G)
end

function exact_psd_factor(P::AbstractMatrix{T}, ::Type{R}) where {T,R}
  work = Matrix{T}(P)
  n = size(work, 1)
  columns = Vector{T}[]
  weights = R[]
  for _ in 1:n
    diagonal = R[real(work[i, i]) for i in 1:n]
    any(<(0), diagonal) &&
      throw(ArgumentError("the exact canonical dark residual is not positive semidefinite"))
    pivot = findfirst(>(0), diagonal)
    if pivot === nothing
      iszero(work) || throw(
        ArgumentError("the exact canonical dark residual is not positive semidefinite")
      )
      break
    end
    weight = diagonal[pivot]
    column = work[:, pivot] / weight
    push!(columns, column)
    push!(weights, weight)
    work -= weight * column * adjoint(column)
  end
  newborn = isempty(columns) ? exact_zeros(T, n, 0) : reduce(hcat, columns)
  return newborn, weights
end

function exact_tangent_lift(
  B::AbstractMatrix{T},
  weights::AbstractVector{R},
  target::AbstractMatrix{T},
  G::AbstractMatrix{T},
) where {T,R}
  n, k = size(B)
  k == 0 && return exact_zeros(T, n, 0)
  W = LinearAlgebra.Diagonal(Vector{T}(weights))
  Bw = B * W
  equations(Y) = vcat(
    exact_hermitian_coordinates(Bw * adjoint(Y) + Y * adjoint(Bw), R),
    exact_hermitian_coordinates(
      im * (W * adjoint(B) * G * Y * W - adjoint(W * adjoint(B) * G * Y * W)) / 2, R
    ),
  )
  unknowns = 2 * n * k
  A = exact_zeros(R, length(equations(exact_zeros(T, n, k))), unknowns)
  for u in 1:unknowns
    Y = exact_zeros(T, n, k)
    index = (u + 1) ÷ 2
    Y[index] = isodd(u) ? one(T) : im * one(T)
    A[:, u] = equations(Y)
  end
  rhs = vcat(exact_hermitian_coordinates(target, R), exact_zeros(R, k^2))
  x = exact_min_norm_solve(A, rhs, Matrix{R}(LinearAlgebra.I, unknowns, unknowns))
  Y = exact_zeros(T, n, k)
  for u in 1:unknowns
    index = (u + 1) ÷ 2
    Y[index] += isodd(u) ? x[u] : im * x[u]
  end
  return Y
end

function native_exact_static_solve(
  residual::Matrix{T},
  known::Matrix{T},
  active::Matrix{T},
  weights::Vector{R},
  gauge_images::Vector{Matrix{T}},
  metric::Matrix{T},
  gauge_metric::Matrix{R},
) where {T,R}
  n = size(residual, 1)
  size(known) == size(metric) == (n, n) && size(active, 1) == n || throw(
    DimensionMismatch("exact static slot inputs must share the Kossakowski dimension")
  )
  length(weights) == size(active, 2) ||
    throw(DimensionMismatch("one rate weight is required per active channel"))
  size(gauge_metric) == (length(gauge_images), length(gauge_images)) ||
    throw(DimensionMismatch("the gauge metric must match the number of gauge directions"))

  Q = exact_dark_projector(active, metric)
  units = exact_hermitian_units(T, n)
  M = exact_form_metric(metric, units, R)
  dark(X) = exact_hermitian_coordinates(Q * X * adjoint(Q), R)
  b = dark(residual - known)
  A =
    isempty(gauge_images) ? exact_zeros(R, length(b), 0) : reduce(hcat, dark.(gauge_images))
  coordinates = exact_min_norm_solve(
    Matrix(transpose(A)) * M * A, -(Matrix(transpose(A)) * M * b), gauge_metric
  )

  solved = copy(residual)
  for (coordinate, image) in zip(coordinates, gauge_images)
    solved += coordinate * image
  end
  P = Q * (solved - known) * adjoint(Q)
  newborn, newborn_weights = exact_psd_factor(P, R)
  born = newborn * LinearAlgebra.Diagonal(Vector{T}(newborn_weights)) * adjoint(newborn)
  correction = exact_tangent_lift(active, weights, solved - known - born, metric)
  Bw = active * LinearAlgebra.Diagonal(Vector{T}(weights))
  coefficient = known + Bw * adjoint(correction) + correction * adjoint(Bw) + born
  coefficient == solved ||
    throw(ArgumentError("exact native static slot failed to reconstruct its coefficient"))
  return ExactStaticSolution(
    coordinates, newborn, newborn_weights, correction, coefficient, P
  )
end

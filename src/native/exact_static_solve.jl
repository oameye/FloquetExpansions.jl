struct NativePositivityError <: Exception
  message::String
end

Base.showerror(io::IO, error::NativePositivityError) = print(io, error.message)

struct ExactStaticSolution{T,R}
  coordinates::Vector{R}
  newborn::Matrix{T}
  newborn_weights::Vector{R}
  correction::Matrix{T}
  coefficient::Matrix{T}
  dark_residual::Matrix{T}
end

exact_zeros(::Type{S}, dims::Int...) where {S} = fill!(Array{S}(undef, dims...), zero(S))

function exact_columns(::Type{S}, height::Int, columns::AbstractVector) where {S}
  result = exact_zeros(S, height, length(columns))
  for (j, column) in pairs(columns)
    result[:, j] = column
  end
  return result
end

function exact_find_pivot(M::AbstractMatrix, from::Int, column::Int)
  for r in from:size(M, 1)
    iszero(M[r, column]) || return r
  end
  return 0
end

function exact_swap_rows!(M::AbstractMatrix, i::Int, j::Int)
  i == j && return M
  for c in axes(M, 2)
    M[i, c], M[j, c] = M[j, c], M[i, c]
  end
  return M
end

function exact_scale_row!(M::AbstractMatrix, row::Int, scale)
  for c in axes(M, 2)
    M[row, c] *= scale
  end
  return M
end

# M[row, :] -= factor * M[pivot, :]
function exact_row_axpy!(M::AbstractMatrix, row::Int, factor, pivot::Int)
  for c in axes(M, 2)
    M[row, c] -= factor * M[pivot, c]
  end
  return M
end

# Explicit-loop products keep the exact path off the wrapper-typed, size-specialised
# LinearAlgebra matmul kernels, which cost far more to compile than to run on these
# small exact matrices.
function exact_mul(A::AbstractMatrix{S}, B::AbstractMatrix{S}) where {S}
  size(A, 2) == size(B, 1) || throw(DimensionMismatch("exact product dimension mismatch"))
  result = exact_zeros(S, size(A, 1), size(B, 2))
  for j in axes(B, 2), k in axes(A, 2)
    b = B[k, j]
    iszero(b) && continue
    for i in axes(A, 1)
      result[i, j] += A[i, k] * b
    end
  end
  return result
end

function exact_mul(A::AbstractMatrix{S}, x::AbstractVector{S}) where {S}
  size(A, 2) == length(x) || throw(DimensionMismatch("exact product dimension mismatch"))
  result = exact_zeros(S, size(A, 1))
  for k in axes(A, 2)
    b = x[k]
    iszero(b) && continue
    for i in axes(A, 1)
      result[i] += A[i, k] * b
    end
  end
  return result
end

exact_mul(A, B, C, rest...) = exact_mul(exact_mul(A, B), C, rest...)

function exact_adjoint(A::AbstractMatrix{S}) where {S}
  return S[conj(A[i, j]) for j in axes(A, 2), i in axes(A, 1)]
end

function exact_diagonal(::Type{S}, values::AbstractVector) where {S}
  result = exact_zeros(S, length(values), length(values))
  for (i, v) in pairs(values)
    result[i, i] = v
  end
  return result
end

function exact_solve(A::AbstractMatrix{S}, B::AbstractVecOrMat{S}) where {S}
  n = size(A, 1)
  size(A, 2) == n || throw(DimensionMismatch("exact solve needs a square matrix"))
  M = Matrix{S}(A)
  X = Array{S}(B)
  Xm = reshape(X, n, :)
  for column in 1:n
    pivot = exact_find_pivot(M, column, column)
    pivot == 0 && throw(ArgumentError("exact solve met a singular matrix"))
    exact_swap_rows!(M, column, pivot)
    exact_swap_rows!(Xm, column, pivot)
    scale = inv(M[column, column])
    exact_scale_row!(M, column, scale)
    exact_scale_row!(Xm, column, scale)
    for r in 1:n
      factor = M[r, column]
      (r == column || iszero(factor)) && continue
      exact_row_axpy!(M, r, factor, column)
      exact_row_axpy!(Xm, r, factor, column)
    end
  end
  return X
end

function exact_row_basis(A::AbstractMatrix{T}) where {T}
  basis = Vector{T}[]
  pivots = Int[]
  selected = Int[]
  for r in axes(A, 1)
    v = T[A[r, c] for c in axes(A, 2)]
    for (b, p) in zip(basis, pivots)
      iszero(v[p]) && continue
      factor = v[p] / b[p]
      for c in eachindex(v)
        v[c] -= factor * b[c]
      end
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
  solution = exact_mul(weighted, exact_solve(exact_mul(Ar, weighted), b[rows]))
  exact_mul(A, solution) == b || throw(ArgumentError("exact linear system is inconsistent"))
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
  Bd = exact_adjoint(B)
  return Matrix{T}(LinearAlgebra.I, n, n) -
         exact_mul(B, exact_solve(exact_mul(Bd, G, B), exact_mul(Bd, G)))
end

struct ExactDarkMap{T}
  B::Matrix{T}
  C::Matrix{T}
end

function exact_dark_map(B::AbstractMatrix{T}, G::AbstractMatrix{T}) where {T}
  size(B, 2) == 0 && return ExactDarkMap{T}(Matrix{T}(B), exact_zeros(T, 0, size(G, 1)))
  Bd = exact_adjoint(B)
  BdG = exact_mul(Bd, G)
  return ExactDarkMap{T}(Matrix{T}(B), exact_solve(exact_mul(BdG, B), BdG))
end

function dark_sandwich(map::ExactDarkMap{T}, X::AbstractMatrix{T}) where {T}
  size(map.C, 1) == 0 && return Matrix(X)
  Y = X - exact_mul(map.B, exact_mul(map.C, X))
  return Y - exact_mul(Y, exact_adjoint(map.C), exact_adjoint(map.B))
end

function exact_isidentity(G::AbstractMatrix)
  return all(
    G[i, j] == (i == j ? one(eltype(G)) : zero(eltype(G))) for
    j in axes(G, 2), i in axes(G, 1)
  )
end

function exact_nonzeros(X::AbstractMatrix)
  return [(a, b, X[a, b]) for b in axes(X, 2) for a in axes(X, 1) if !iszero(X[a, b])]
end

function exact_trace_product(entries, Y::AbstractMatrix{T}) where {T}
  value = zero(T)
  for (a, b, x) in entries
    y = Y[b, a]
    iszero(y) || (value += x * y)
  end
  return real(value)
end

function exact_psd_factor(P::AbstractMatrix{T}, ::Type{R}) where {T,R}
  work = Matrix{T}(P)
  n = size(work, 1)
  columns = Vector{T}[]
  weights = R[]
  for _ in 1:n
    diagonal = R[real(work[i, i]) for i in 1:n]
    any(<(0), diagonal) && throw(
      NativePositivityError(
        "the exact canonical dark residual is not positive semidefinite"
      ),
    )
    pivot = findfirst(>(0), diagonal)
    if pivot === nothing
      iszero(work) || throw(
        NativePositivityError(
          "the exact canonical dark residual is not positive semidefinite"
        ),
      )
      break
    end
    weight = diagonal[pivot]
    column = T[work[i, pivot] / weight for i in 1:n]
    push!(columns, column)
    push!(weights, weight)
    for j in 1:n
      scaled = weight * conj(column[j])
      iszero(scaled) && continue
      for i in 1:n
        work[i, j] -= column[i] * scaled
      end
    end
  end
  newborn = exact_columns(T, n, columns)
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
  Bd = exact_adjoint(B)
  BdG = exact_mul(Bd, G)
  H = exact_mul(BdG, B)
  C = exact_solve(H, BdG)
  W = exact_diagonal(T, weights)
  Winv = exact_diagonal(T, inv.(weights))
  Ct = exact_adjoint(C)
  dark = exact_mul(target - exact_mul(B, C, target), Ct, Winv)
  M = exact_mul(C, target, Ct)
  system = exact_zeros(R, 4 * k^2, 2 * k^2)
  for u in 1:(2 * k ^ 2)
    A = exact_zeros(T, k, k)
    A[(u + 1) ÷ 2] = isodd(u) ? one(T) : im * one(T)
    system[:, u] = exact_tangent_equations(A, W, H, R)
  end
  unknowns = 2 * k^2
  rhs = vcat(exact_complex_coordinates(M, R), exact_zeros(R, 2 * k^2))
  x = exact_min_norm_solve(system, rhs, Matrix{R}(LinearAlgebra.I, unknowns, unknowns))
  A = exact_zeros(T, k, k)
  for u in 1:unknowns
    A[(u + 1) ÷ 2] += isodd(u) ? x[u] : im * x[u]
  end
  return exact_mul(B, A) + dark
end

function exact_tangent_equations(
  A::Matrix{T}, W::Matrix{T}, H::Matrix{T}, ::Type{R}
) where {T,R}
  WHAW = exact_mul(W, H, A, W)
  return vcat(
    exact_complex_coordinates(exact_mul(W, exact_adjoint(A)) + exact_mul(A, W), R),
    exact_complex_coordinates(WHAW - exact_adjoint(WHAW), R),
  )
end

function exact_complex_coordinates(X::AbstractMatrix, ::Type{R}) where {R}
  return vcat(R[real(x) for x in vec(X)], R[imag(x) for x in vec(X)])
end

function exact_combination(J::Vector{Matrix{T}}, x::AbstractVector, n::Int) where {T}
  result = exact_zeros(T, n, n)
  for (k, xk) in pairs(x)
    iszero(xk) && continue
    result += xk * J[k]
  end
  return result
end

function exact_column_product(F::AbstractMatrix{T}, columns, N::AbstractMatrix) where {T}
  result = exact_zeros(T, size(F, 1), size(N, 2))
  for (j, column) in pairs(columns), a in axes(F, 1)
    f = F[a, column]
    iszero(f) && continue
    for l in axes(N, 2)
      v = N[j, l]
      iszero(v) || (result[a, l] += f * v)
    end
  end
  return result
end

function exact_row_cancellation(P0, J, C, c, ::Type{R}) where {R}
  n = size(P0, 1)
  for i in 1:n
    rows = exact_columns(R, 2 * n, [exact_complex_coordinates(Jk[:, i:i], R) for Jk in J])
    rhs = -exact_complex_coordinates(P0[:, i:i], R)
    iszero(rows) && iszero(rhs) && continue
    trialC = vcat(C, rows)
    trialc = vcat(c, rhs)
    exact_consistent_solve(trialC, trialc, size(C, 2)) === nothing ||
      ((C, c) = (trialC, trialc))
  end
  return C, c
end

function exact_fixed_block(free, n::Int)
  fixed = [i for i in 1:n if all(iszero(F[i, i]) for F in free)]
  while !isempty(fixed)
    violations = [count(j -> any(!iszero(F[i, j]) for F in free), fixed) for i in fixed]
    worst = argmax(violations)
    violations[worst] == 0 && return fixed
    deleteat!(fixed, worst)
  end
  return fixed
end

function exact_free_directions(J::Vector{Matrix{T}}, Z::Matrix, n::Int) where {T}
  Z == LinearAlgebra.I && return J
  free = Matrix{T}[]
  for j in axes(Z, 2)
    push!(free, exact_combination(J, view(Z, :, j), n))
  end
  return free
end

function exact_coordinate_columns(
  matrices::Vector{Matrix{T}}, fixed::Vector{Int}, N::Matrix{T}, ::Type{R}, height::Int
) where {T,R}
  result = exact_zeros(R, height, length(matrices))
  for (j, F) in pairs(matrices)
    result[:, j] = exact_complex_coordinates(exact_column_product(F, fixed, N), R)
  end
  return result
end

function exact_facial_step(
  P0::Matrix{T}, J::Vector{Matrix{T}}, x0::Vector{R}, Z::Matrix{R}, ::Type{R}
) where {T,R}
  n = size(P0, 1)
  free = exact_free_directions(J, Z, n)
  fixed = exact_fixed_block(free, n)
  isempty(fixed) && return nothing
  current = P0 + exact_combination(J, x0, n)
  block = current[fixed, fixed]
  exact_is_psd(block, R) ||
    throw(NativePositivityError("no PSD lift exists in this static gauge family"))
  kernel = exact_nullspace(block)
  isempty(kernel) && return nothing
  N = exact_columns(T, length(fixed), kernel)
  height = 2 * n * size(N, 2)
  rows = exact_coordinate_columns(free, fixed, N, R, height)
  rhs = -exact_complex_coordinates(exact_mul(current[:, fixed], N), R)
  y = exact_consistent_solve(rows, rhs, size(Z, 2))
  y === nothing &&
    throw(NativePositivityError("no PSD lift exists in this static gauge family"))
  constraint = free === J ? rows : exact_coordinate_columns(J, fixed, N, R, height)
  constraint_rhs = -exact_complex_coordinates(exact_mul(P0[:, fixed], N), R)
  return constraint, constraint_rhs, y, exact_nullspace(rows)
end

function exact_consistent_solve(A::Matrix{R}, b::Vector{R}, unknowns::Int) where {R}
  try
    return exact_min_norm_solve(A, b, Matrix{R}(LinearAlgebra.I, unknowns, unknowns))
  catch error
    error isa ArgumentError || rethrow()
    return nothing
  end
end

function exact_facial_constraints(
  P0::Matrix{T}, J::Vector{Matrix{T}}, ::Type{R}
) where {T,R}
  p = length(J)
  x0 = exact_zeros(R, p)
  Z = Matrix{R}(LinearAlgebra.I, p, p)
  C = exact_zeros(R, 0, p)
  c = R[]
  while true
    step = exact_facial_step(P0, J, x0, Z, R)
    step === nothing && break
    constraint, rhs, y, reduced = step
    C = vcat(C, constraint)
    c = vcat(c, rhs)
    x0 += exact_mul(Z, y)
    length(reduced) == size(Z, 2) && break
    Z = exact_mul(Z, exact_columns(R, size(Z, 2), reduced))
  end
  return C, c, Z
end

function exact_section(C, c, Z, normal_matrix, normal_rhs, gauge_metric)
  Zt = Matrix(transpose(Z))
  return exact_min_norm_solve(
    vcat(C, exact_mul(Zt, normal_matrix)), vcat(c, exact_mul(Zt, normal_rhs)), gauge_metric
  )
end

function exact_is_psd(P::AbstractMatrix, ::Type{R}) where {R}
  try
    exact_psd_factor(P, R)
    return true
  catch error
    error isa NativePositivityError || rethrow()
    return false
  end
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

  dark = exact_dark_map(active, metric)
  P0 = dark_sandwich(dark, residual - known)
  J = [dark_sandwich(dark, image) for image in gauge_images]
  identity_metric = exact_isidentity(metric)::Bool
  GJ = identity_metric ? J : [exact_mul(metric, X) for X in J]
  GP = identity_metric ? P0 : exact_mul(metric, P0)
  entries = [exact_nonzeros(X) for X in GJ]
  p = length(J)
  normal_matrix = exact_zeros(R, p, p)
  for i in 1:p, j in i:p
    value = exact_trace_product(entries[i], GJ[j])
    normal_matrix[i, j] = value
    normal_matrix[j, i] = value
  end
  normal_rhs = R[-exact_trace_product(entries[i], GP) for i in 1:p]
  C, c, Z = exact_facial_constraints(P0, J, R)
  coordinates = exact_section(C, c, Z, normal_matrix, normal_rhs, gauge_metric)
  if !exact_is_psd(P0 + exact_combination(J, coordinates, n), R)
    C, c = exact_row_cancellation(P0, J, C, c, R)
    kernel = exact_nullspace(C)
    Z = exact_columns(R, p, kernel)
    coordinates = exact_section(C, c, Z, normal_matrix, normal_rhs, gauge_metric)
  end

  solved = copy(residual)
  for (coordinate, image) in zip(coordinates, gauge_images)
    solved += coordinate * image
  end
  P = dark_sandwich(dark, solved - known)
  newborn, newborn_weights = exact_psd_factor(P, R)
  born = exact_mul(newborn, exact_diagonal(T, newborn_weights), exact_adjoint(newborn))
  correction = exact_tangent_lift(active, weights, solved - known - born, metric)
  Bw = exact_mul(active, exact_diagonal(T, weights))
  coefficient =
    known +
    exact_mul(Bw, exact_adjoint(correction)) +
    exact_mul(correction, exact_adjoint(Bw)) +
    born
  coefficient == solved ||
    throw(ArgumentError("exact native static slot failed to reconstruct its coefficient"))
  return ExactStaticSolution(
    coordinates, newborn, newborn_weights, correction, coefficient, P
  )
end

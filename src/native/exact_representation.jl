struct ExactLiouvilleRepresentation{T,R} <: NativeRepresentation
  d::Int
  hilbert::Vector{R}
  frame::Vector{Matrix{T}}
  full::Vector{Matrix{T}}
  inverse_frame::Matrix{T}
  metric::Matrix{T}
  gauge::Vector{Matrix{T}}
  gauge_metric::Matrix{R}
end

function ExactLiouvilleRepresentation{T}(d::Int) where {T}
  return ExactLiouvilleRepresentation{T}(ones(real(T), d))
end

function ExactLiouvilleRepresentation{T}(hilbert::Vector{R}) where {T,R<:Real}
  d = length(hilbert)
  d >= 2 || throw(ArgumentError("an exact Liouville representation needs dimension d >= 2"))
  all(>(0), hilbert) ||
    throw(ArgumentError("the Hilbert metric of an exact representation must be positive"))
  g = Vector{real(T)}(hilbert)
  frame = exact_traceless_frame(T, d)
  full = vcat([exact_identity(T, d)], frame)
  V = reduce(hcat, [vec(F) for F in full])
  inverse_frame = exact_solve(V, Matrix{T}(LinearAlgebra.I, d^2, d^2))
  metric = T[LinearAlgebra.tr(hilbert_adjoint(g, F) * G) for F in frame, G in frame]
  gauge = exact_gauge_generators(g, frame)
  gauge_metric = superoperator_metric(g, gauge)
  return ExactLiouvilleRepresentation{T,real(T)}(
    d, g, frame, full, inverse_frame, metric, gauge, gauge_metric
  )
end

function hilbert_adjoint(g::AbstractVector, X::AbstractMatrix{T}) where {T}
  d = length(g)
  return T[conj(X[b, a]) * g[b] / g[a] for a in 1:d, b in 1:d]
end

function superoperator_metric(g::AbstractVector{R}, directions) where {R}
  d = length(g)
  weight = R[g[(p - 1) % d + 1] / g[(p - 1) ÷ d + 1] for p in 1:(d ^ 2)]
  supports = [findall(!iszero, A) for A in directions]
  n = length(directions)
  metric = exact_zeros(R, n, n)
  for k in 1:n, l in k:n
    value = zero(R)
    A = directions[k]
    B = directions[l]
    for index in supports[k]
      q, p = Tuple(index)
      iszero(B[q, p]) && continue
      value += real(conj(A[q, p]) * B[q, p]) * weight[q] / weight[p]
    end
    metric[k, l] = value
    metric[l, k] = value
  end
  return metric
end

exact_identity(::Type{T}, d::Int) where {T} = Matrix{T}(LinearAlgebra.I, d, d)
exact_left(A::AbstractMatrix{T}) where {T} = kron(exact_identity(T, size(A, 1)), A)
function exact_right(B::AbstractMatrix{T}) where {T}
  return kron(transpose(B), exact_identity(T, size(B, 1)))
end
function exact_sandwich(g::AbstractVector, X::AbstractMatrix, Y::AbstractMatrix)
  return kron(transpose(hilbert_adjoint(g, Y)), X)
end

function exact_unit(::Type{T}, d::Int, a::Int, b::Int) where {T}
  E = exact_zeros(T, d, d)
  E[a, b] = one(T)
  return E
end

function exact_traceless_frame(::Type{T}, d::Int) where {T}
  frame = Matrix{T}[]
  for b in 1:d, a in 1:d
    a == b || push!(frame, exact_unit(T, d, a, b))
  end
  for a in 1:(d - 1)
    push!(frame, exact_unit(T, d, a, a) - exact_identity(T, d) / d)
  end
  return frame
end

function exact_hermitian_frame(g::AbstractVector, frame::Vector{Matrix{T}}) where {T}
  hermitian = Matrix{T}[]
  d = length(g)
  for F in frame
    Fd = hilbert_adjoint(g, F)
    if Fd == F
      push!(hermitian, F)
    else
      a, b = Tuple(findfirst(!iszero, F))
      a < b && push!(hermitian, F + Fd, im * (F - Fd))
    end
  end
  return hermitian
end

function exact_dissipator(
  g::AbstractVector, frame::Vector{Matrix{T}}, C::AbstractMatrix{T}
) where {T}
  n = size(frame[1], 1)^2
  result = exact_zeros(T, n, n)
  for j in eachindex(frame), i in eachindex(frame)
    iszero(C[i, j]) && continue
    F = hilbert_adjoint(g, frame[j]) * frame[i]
    result +=
      C[i, j] *
      (exact_sandwich(g, frame[i], frame[j]) - exact_left(F) / 2 - exact_right(F) / 2)
  end
  return result
end

function exact_gauge_generators(g::AbstractVector, frame::Vector{Matrix{T}}) where {T}
  generators = Matrix{T}[]
  for h in exact_hermitian_frame(g, frame)
    push!(generators, -im * (exact_left(h) - exact_right(h)))
  end
  for unit in exact_hermitian_units(T, length(frame))
    push!(generators, exact_dissipator(g, frame, unit))
  end
  return generators
end

function exact_chi_matrix(representation::ExactLiouvilleRepresentation{T}, L) where {T}
  d = representation.d
  g = representation.hilbert
  units = exact_zeros(T, d^2, d^2)
  for f in 1:d, e in 1:d, b in 1:d, a in 1:d
    units[(b - 1) * d + a, (f - 1) * d + e] =
      L[(e - 1) * d + a, (f - 1) * d + b] * g[f] / g[e]
  end
  W = representation.inverse_frame
  return W * units * adjoint(W)
end

function native_kossakowski(representation::ExactLiouvilleRepresentation, L::AbstractMatrix)
  return exact_chi_matrix(representation, L)[2:end, 2:end]
end

function native_hamiltonian(representation::ExactLiouvilleRepresentation{T}, L) where {T}
  chi = exact_chi_matrix(representation, L)
  G = exact_zeros(T, representation.d, representation.d)
  for (i, F) in pairs(representation.frame)
    G += chi[i + 1, 1] * F
  end
  return (im // 2) * (G - hilbert_adjoint(representation.hilbert, G))
end

function native_gksl(
  representation::ExactLiouvilleRepresentation{T}, H::AbstractMatrix, C::AbstractMatrix
) where {T}
  Hm = Matrix{T}(H)
  return -im * (exact_left(Hm) - exact_right(Hm)) +
         exact_dissipator(representation.hilbert, representation.frame, Matrix{T}(C))
end

native_gauge_directions(representation::ExactLiouvilleRepresentation) = representation.gauge

native_matches(::ExactLiouvilleRepresentation, A, B, tol::Real) = A == B

function native_dark_target(
  representation::ExactLiouvilleRepresentation, residual, known, active, tol::Real
)
  Q = exact_dark_projector(active, representation.metric)
  return Q * (residual - known) * adjoint(Q)
end

function native_slot_solve(
  representation::ExactLiouvilleRepresentation{T,R},
  residual,
  known,
  active,
  weights,
  images,
  directions,
  tol::Real,
) where {T,R}
  gauge_metric = superoperator_metric(representation.hilbert, directions)
  solution = native_exact_static_solve(
    residual,
    known,
    active,
    R[real(weight) for weight in weights],
    images,
    representation.metric,
    gauge_metric,
  )
  return solution, solution.correction, T.(solution.newborn_weights)
end

struct NativeFrameError <: Exception
  message::String
end

Base.showerror(io::IO, error::NativeFrameError) = print(io, error.message)

struct AlgebraicLiouvilleRepresentation{T,R} <: NativeRepresentation
  algebra::OperatorAlgebra{T,R}
  frame::Vector{Monomial}
  index::Dict{Monomial,Int}
  right_basis::Matrix{T}
  right_identity::Vector{T}
  metric::Matrix{T}
  gauge_degree::Int
end

function site_monomials(site::AlgebraSite, degree::Int)
  if site.kind == LEVEL_SITE
    N = site.levels
    monomials = [(0, 0)]
    for i in 1:N, j in 1:N
      (i, j) == (N, N) || push!(monomials, (i, j))
    end
    return monomials
  end
  if site.kind == SPIN_SITE
    return [
      (a, b, total - a - b) for total in 0:degree for a in 0:total for b in 0:(total - a)
    ]
  end
  return [(k, total - k) for total in 0:degree for k in 0:total]
end

function algebra_monomials(algebra::OperatorAlgebra, degree::Int)
  monomials = Tuple{Monomial,Int}[(Int[], 0)]
  for site in algebra.sites
    next = Tuple{Monomial,Int}[]
    for (monomial, partial) in monomials, slots in site_monomials(site, degree)
      total = site.kind == LEVEL_SITE ? partial : partial + sum(slots)
      total <= degree && push!(next, (append!(copy(monomial), slots), total))
    end
    monomials = next
  end
  identity = identity_monomial(algebra)
  return sort!(Monomial[m for (m, _) in monomials if m != identity])
end

function AlgebraicLiouvilleRepresentation(
  algebra::OperatorAlgebra{T,R}, degree::Int, gauge_degree::Int
) where {T,R}
  frame = algebra_monomials(algebra, degree)
  index = Dict(m => k for (k, m) in pairs(frame))
  n = length(frame)
  identity = identity_monomial(algebra)
  basis = exact_zeros(T, n + 1, n + 1)
  basis[n + 1, n + 1] = one(T)
  for (ν, F) in pairs(frame)
    for (M, c) in monomial_adjoint(algebra, F)
      row = M == identity ? n + 1 : get(index, M, 0)
      row == 0 &&
        throw(ArgumentError("the algebraic frame is not closed under the adjoint"))
      basis[row, ν] += c
    end
  end
  coordinates = exact_solve(basis, Matrix{T}(LinearAlgebra.I, n + 1, n + 1))
  right_basis = coordinates[1:n, 1:n]
  right_identity = coordinates[n + 1, 1:n]
  metric = Matrix{T}(LinearAlgebra.I, n, n)
  return AlgebraicLiouvilleRepresentation{T,R}(
    algebra, frame, index, right_basis, right_identity, metric, gauge_degree
  )
end

function frame_operator(representation::AlgebraicLiouvilleRepresentation, μ::Int)
  return algebra_operator(representation.algebra, [(representation.frame[μ], 1)])
end

function frame_coordinates(
  representation::AlgebraicLiouvilleRepresentation{T}, X::AlgebraOperator
) where {T}
  coordinates = exact_zeros(T, length(representation.frame))
  identity = identity_monomial(representation.algebra)
  for (m, c) in X.terms
    m == identity && continue
    k = get(representation.index, m, 0)
    k == 0 && throw(NativeFrameError("operator has a monomial outside the algebraic frame"))
    coordinates[k] = c
  end
  return coordinates
end

function frame_operator(
  representation::AlgebraicLiouvilleRepresentation, coordinates::AbstractVector
)
  return algebra_operator(
    representation.algebra,
    [(representation.frame[k], c) for (k, c) in pairs(coordinates) if !iszero(c)],
  )
end

function native_kossakowski(
  representation::AlgebraicLiouvilleRepresentation{T}, L::AlgebraSuperoperator
) where {T}
  n = length(representation.frame)
  identity = identity_monomial(representation.algebra)
  C = exact_zeros(T, n, n)
  for ((X, M), c) in L.terms
    (X == identity || M == identity) && continue
    μ = get(representation.index, X, 0)
    k = get(representation.index, M, 0)
    (μ == 0 || k == 0) && throw(
      NativeFrameError("superoperator has a sandwich term outside the algebraic frame")
    )
    for ν in 1:n
      w = representation.right_basis[ν, k]
      iszero(w) || (C[μ, ν] += c * w)
    end
  end
  return C
end

function native_hamiltonian(
  representation::AlgebraicLiouvilleRepresentation{T,R}, L::AlgebraSuperoperator
) where {T,R}
  algebra = representation.algebra
  identity = identity_monomial(algebra)
  K = Dict{Monomial,T}()
  for ((X, M), c) in L.terms
    if M == identity
      accumulate!(K, X, c)
    elseif X != identity
      k = get(representation.index, M, 0)
      k == 0 &&
        throw(NativeFrameError("superoperator has a term outside the algebraic frame"))
      accumulate!(K, X, c * representation.right_identity[k])
    end
  end
  G = AlgebraOperator{T,R}(algebra, K)
  H = (im // 2) * (G - adjoint(G))
  delete!(H.terms, identity)
  return H
end

function native_gksl(
  representation::AlgebraicLiouvilleRepresentation{T,R},
  H::AlgebraOperator,
  C::AbstractMatrix,
) where {T,R}
  result = -im * (left_action(H) - right_action(H))
  n = length(representation.frame)
  for ν in 1:n, μ in 1:n
    iszero(C[μ, ν]) && continue
    F = frame_operator(representation, μ)
    Gd = adjoint(frame_operator(representation, ν))
    norm = Gd * F
    result =
      result +
      T(C[μ, ν]) *
      (sandwich(F, Gd) - (1 // 2) * left_action(norm) - (1 // 2) * right_action(norm))
  end
  return result
end

native_matches(::AlgebraicLiouvilleRepresentation, A, B, ::Real) = A == B

function native_dark_target(
  representation::AlgebraicLiouvilleRepresentation, residual, known, active, ::Real
)
  return dark_sandwich(exact_dark_map(active, representation.metric), residual - known)
end

function native_slot_solve(
  representation::AlgebraicLiouvilleRepresentation{T,R},
  residual,
  known,
  active,
  weights,
  family,
  ::Real,
) where {T,R}
  solution = native_exact_static_solve(
    residual,
    known,
    active,
    R[real(weight) for weight in weights],
    family.images,
    representation.metric,
    family.metric,
  )
  return solution, solution.correction, T.(solution.newborn_weights)
end

function native_gauge_metric(
  ::AlgebraicLiouvilleRepresentation{T,R},
  directions::AbstractVector{<:AlgebraSuperoperator},
) where {T,R}
  n = length(directions)
  metric = exact_zeros(R, n, n)
  for k in 1:n, l in k:n
    A = directions[k].terms
    B = directions[l].terms
    value = zero(T)
    for (key, v) in B
      a = get(A, key, zero(T))
      iszero(a) || (value += conj(a) * v)
    end
    metric[k, l] = real(value)
    metric[l, k] = real(value)
  end
  return metric
end

function hermitian_generators(
  representation::AlgebraicLiouvilleRepresentation{T,R}, @nospecialize(keep)
) where {T,R}
  generators = AlgebraOperator{T,R}[]
  seen = Set{Monomial}()
  for (μ, F) in pairs(representation.frame)
    (F in seen || !(keep(μ, μ)::Bool)) && continue
    X = frame_operator(representation, μ)
    Xd = adjoint(X)
    for (m, _) in Xd.terms
      push!(seen, m)
    end
    push!(generators, X + Xd)
    isequal(X, Xd) || push!(generators, im * (X - Xd))
  end
  return generators
end

# Restricts a direction filter to frame indices of bounded degree. A callable struct rather
# than a closure, so the captured filter stays a concrete type parameter.
struct RestrictedKeep{K}
  small::Vector{Int}
  keep::K
end

(r::RestrictedKeep)(μ::Int, ν::Int) = μ in r.small && r.keep(μ, ν)::Bool

function algebraic_gauge_directions(
  representation::AlgebraicLiouvilleRepresentation{T,R}, @nospecialize(keep)
) where {T,R}
  return algebraic_gauge_directions(representation, keep, representation.gauge_degree)
end

function algebraic_gauge_directions(
  representation::AlgebraicLiouvilleRepresentation{T,R},
  @nospecialize(keep),
  gauge_degree::Int,
) where {T,R}
  algebra = representation.algebra
  n = length(representation.frame)
  degree(μ) = monomial_degree(algebra, representation.frame[μ])
  small = Int[μ for μ in 1:n if degree(μ) <= gauge_degree]
  directions = AlgebraSuperoperator{T,R}[]
  zero_frame = exact_zeros(T, n, n)
  for h in hermitian_generators(representation, RestrictedKeep(small, keep))
    push!(directions, native_gksl(representation, h, zero_frame))
  end
  empty_hamiltonian = AlgebraOperator{T,R}(algebra, Dict{Monomial,T}())
  for (i, μ) in pairs(small), ν in small[i:end]
    keep(μ, ν)::Bool || continue
    if μ == ν
      unit = copy(zero_frame)
      unit[μ, μ] = one(T)
      push!(directions, native_gksl(representation, empty_hamiltonian, unit))
    else
      real_unit = copy(zero_frame)
      real_unit[μ, ν] = one(T)
      real_unit[ν, μ] = one(T)
      imag_unit = copy(zero_frame)
      imag_unit[μ, ν] = im
      imag_unit[ν, μ] = -im
      push!(directions, native_gksl(representation, empty_hamiltonian, real_unit))
      push!(directions, native_gksl(representation, empty_hamiltonian, imag_unit))
    end
  end
  return directions
end

function native_gauge_directions(representation::AlgebraicLiouvilleRepresentation)
  return algebraic_gauge_directions(representation, (_, _) -> true)
end

function native_gauge_families(
  representation::AlgebraicLiouvilleRepresentation, ::NoHomologicalInverse, L0
)
  return [
    native_gauge_family(
      representation, algebraic_gauge_directions(representation, (_, _) -> true, degree), L0
    ) for degree in 1:(representation.gauge_degree)
  ]
end

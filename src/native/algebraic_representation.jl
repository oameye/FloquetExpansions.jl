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
  return [(k, total - k) for total in 0:degree for k in 0:total]
end

function algebra_monomials(algebra::OperatorAlgebra, degree::Int)
  monomials = [Int[]]
  for site in algebra.sites
    next = Monomial[]
    for monomial in monomials, (p, q) in site_monomials(site, degree)
      candidate = vcat(monomial, [p, q])
      partial = 0
      for (s, other) in pairs(algebra.sites[1:(length(candidate) ÷ 2)])
        other.kind == LEVEL_SITE || (partial += candidate[2s - 1] + candidate[2s])
      end
      partial <= degree && push!(next, candidate)
    end
    monomials = next
  end
  identity = identity_monomial(algebra)
  return sort!(filter(!=(identity), monomials))
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
    k == 0 && throw(ArgumentError("operator has a monomial outside the algebraic frame"))
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
    (μ == 0 || k == 0) &&
      throw(ArgumentError("superoperator has a sandwich term outside the algebraic frame"))
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
      k = representation.index[M]
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

native_matches(::AlgebraicLiouvilleRepresentation, A, B, tol::Real) = A == B

function native_dark_target(
  representation::AlgebraicLiouvilleRepresentation, residual, known, active, tol::Real
)
  Q = exact_dark_projector(active, representation.metric)
  return Q * (residual - known) * adjoint(Q)
end

function native_slot_solve(
  representation::AlgebraicLiouvilleRepresentation{T,R},
  residual,
  known,
  active,
  weights,
  images,
  directions,
  tol::Real,
) where {T,R}
  gauge_metric = R[
    real(sum(conj(get(A.terms, k, zero(T))) * v for (k, v) in B.terms; init=zero(T))) for
    A in directions, B in directions
  ]
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

function hermitian_generators(representation::AlgebraicLiouvilleRepresentation, keep)
  generators = AlgebraOperator[]
  seen = Set{Monomial}()
  for (μ, F) in pairs(representation.frame)
    (F in seen || !keep(μ, μ)) && continue
    X = frame_operator(representation, μ)
    Xd = adjoint(X)
    for (m, _) in Xd.terms
      push!(seen, m)
    end
    push!(generators, X + Xd)
    X == Xd || push!(generators, im * (X - Xd))
  end
  return generators
end

function algebraic_gauge_directions(
  representation::AlgebraicLiouvilleRepresentation{T,R}, keep
) where {T,R}
  algebra = representation.algebra
  n = length(representation.frame)
  degree(μ) = monomial_degree(algebra, representation.frame[μ])
  small = [μ for μ in 1:n if degree(μ) <= representation.gauge_degree]
  directions = AlgebraSuperoperator{T,R}[]
  zero_frame = exact_zeros(T, n, n)
  for h in hermitian_generators(representation, (μ, ν) -> μ in small && keep(μ, ν))
    push!(directions, native_gksl(representation, h, zero_frame))
  end
  empty_hamiltonian = AlgebraOperator{T,R}(algebra, Dict{Monomial,T}())
  for (i, μ) in pairs(small), ν in small[i:end]
    keep(μ, ν) || continue
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
  return algebraic_gauge_directions(representation, (μ, ν) -> true)
end

function sparse_exact_solve(columns::Vector{Dict{K,T}}, target::Dict{K,T}) where {K,T}
  basis = Dict{K,T}[]
  pivots = K[]
  combinations = Dict{Int,T}[]
  for (j, column) in pairs(columns)
    v = copy(column)
    combination = Dict{Int,T}(j => one(T))
    for (b, p, c) in zip(basis, pivots, combinations)
      factor = get(v, p, zero(T))
      iszero(factor) && continue
      for (k, x) in b
        accumulate!(v, k, -factor * x)
      end
      for (k, x) in c
        accumulate!(combination, k, -factor * x)
      end
    end
    isempty(v) && continue
    p = first(sort!(collect(keys(v))))
    scale = inv(v[p])
    for k in keys(v)
      v[k] *= scale
    end
    for k in keys(combination)
      combination[k] *= scale
    end
    push!(basis, v)
    push!(pivots, p)
    push!(combinations, combination)
  end
  residual = copy(target)
  solution = Dict{Int,T}()
  for (b, p, c) in zip(basis, pivots, combinations)
    factor = get(residual, p, zero(T))
    iszero(factor) && continue
    for (k, x) in b
      accumulate!(residual, k, -factor * x)
    end
    for (k, x) in c
      accumulate!(solution, k, factor * x)
    end
  end
  isempty(residual) || return nothing
  return solution
end

function charge_variables(algebra::OperatorAlgebra)
  offsets = Int[]
  count = 0
  for site in algebra.sites
    push!(offsets, count)
    site.kind == BOSON_SITE && (count += 1)
    site.kind == LEVEL_SITE && (count += site.levels)
  end
  return offsets, count
end

function monomial_charge(
  algebra::OperatorAlgebra, offsets::Vector{Int}, count::Int, M::Monomial
)
  charge = zeros(Int, count)
  for (s, site) in pairs(algebra.sites)
    p, q = M[2s - 1], M[2s]
    if site.kind == BOSON_SITE
      charge[offsets[s] + 1] += p - q
    elseif site.kind == LEVEL_SITE && (p, q) != (0, 0)
      charge[offsets[s] + p] += 1
      charge[offsets[s] + q] -= 1
    end
  end
  return charge
end

struct AlgebraicChargeInverse{T,R} <: HomologicalInverse
  L0::AlgebraSuperoperator{T,R}
  offsets::Vector{Int}
  weights::Matrix{Int}
  directions::Vector{AlgebraSuperoperator{T,R}}
  max_degree::Int
end

function key_label(inverse::AlgebraicChargeInverse, key::Tuple{Monomial,Monomial})
  algebra = inverse.L0.algebra
  count = size(inverse.weights, 2)
  c =
    monomial_charge(algebra, inverse.offsets, count, key[1]) +
    monomial_charge(algebra, inverse.offsets, count, key[2])
  return inverse.weights * c
end

function AlgebraicChargeInverse(
  representation::AlgebraicLiouvilleRepresentation{T,R},
  L0::AlgebraSuperoperator{T,R},
  leading::AbstractMatrix,
  max_degree::Int,
) where {T,R}
  algebra = representation.algebra
  offsets, count = charge_variables(algebra)
  rows = Set{Vector{Int}}()
  for (X, Y) in keys(L0.terms)
    c =
      monomial_charge(algebra, offsets, count, X) +
      monomial_charge(algebra, offsets, count, Y)
    iszero(c) || push!(rows, c)
  end
  for j in axes(leading, 2)
    support = [k for k in axes(leading, 1) if !iszero(leading[k, j])]
    for k in support[2:end]
      c =
        monomial_charge(algebra, offsets, count, representation.frame[k]) -
        monomial_charge(algebra, offsets, count, representation.frame[support[1]])
      iszero(c) || push!(rows, c)
    end
  end
  for (s, site) in pairs(algebra.sites)
    site.kind == LEVEL_SITE || continue
    c = zeros(Int, count)
    c[offsets[s] + 1] = 1
    push!(rows, c)
  end
  kernel = if count == 0
    Vector{Rational{BigInt}}[]
  elseif isempty(rows)
    [Rational{BigInt}.(Matrix{Int}(LinearAlgebra.I, count, count)[:, k]) for k in 1:count]
  else
    exact_nullspace(Rational{BigInt}[row[j] for row in collect(rows), j in 1:count])
  end
  weights = zeros(Int, length(kernel), count)
  for (r, v) in pairs(kernel)
    weights[r, :] = integer_vector(v)
  end
  label(μ) = weights * monomial_charge(algebra, offsets, count, representation.frame[μ])
  directions = algebraic_gauge_directions(representation, (μ, ν) -> label(μ) == label(ν))
  return AlgebraicChargeInverse{T,R}(L0, offsets, weights, directions, max_degree)
end

function charged_ansatz(inverse::AlgebraicChargeInverse, label::Vector{Int}, degree::Int)
  algebra = inverse.L0.algebra
  monomials = vcat([identity_monomial(algebra)], algebra_monomials(algebra, degree))
  keys = Tuple{Monomial,Monomial}[]
  for X in monomials, Y in monomials
    monomial_degree(algebra, X) + monomial_degree(algebra, Y) <= degree || continue
    key_label(inverse, (X, Y)) == label && push!(keys, (X, Y))
  end
  return keys
end

function regular_gauge(
  inverse::AlgebraicChargeInverse{T,R},
  representation::AlgebraicLiouvilleRepresentation{T,R},
  L0,
  residual,
  known,
  active,
  tol,
) where {T,R}
  L0 == inverse.L0 ||
    throw(ArgumentError("charge inverse was built for a different averaged generator"))
  algebra = representation.algebra
  S = zero(L0)
  isempty(inverse.weights) && return S
  target = native_dark_target(representation, residual, known, active, tol)
  empty_hamiltonian = AlgebraOperator{T,R}(algebra, Dict{Monomial,T}())
  Y = native_gksl(representation, empty_hamiltonian, target)
  sectors = Dict{Vector{Int},Dict{Tuple{Monomial,Monomial},T}}()
  for (key, value) in Y.terms
    label = key_label(inverse, key)
    iszero(label) && continue
    accumulate!(get!(sectors, label, Dict{Tuple{Monomial,Monomial},T}()), key, -value)
  end
  for (label, rhs) in sectors
    top = maximum(
      monomial_degree(algebra, X) + monomial_degree(algebra, M) for (X, M) in keys(rhs)
    )
    solution = nothing
    ansatz = Tuple{Monomial,Monomial}[]
    for degree in max(0, top - 1):inverse.max_degree
      ansatz = charged_ansatz(inverse, label, degree)
      columns = [
        (
          L0 * AlgebraSuperoperator{T,R}(algebra, Dict(k => one(T))) -
          AlgebraSuperoperator{T,R}(algebra, Dict(k => one(T))) * L0
        ).terms for k in ansatz
      ]
      solution = sparse_exact_solve(columns, rhs)
      solution === nothing || break
    end
    solution === nothing && throw(
      ArgumentError(
        "charged homological equation has no polynomial solution up to degree $(inverse.max_degree)",
      ),
    )
    for (k, x) in solution
      accumulate!(S.terms, ansatz[k], x)
    end
  end
  return S
end

function singular_gauge_directions(inverse::AlgebraicChargeInverse, ::NativeRepresentation)
  return inverse.directions
end

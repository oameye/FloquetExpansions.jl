const LatticeScalar = Complex{Rational{BigInt}}

# Raised inside the virtual step when the model leaves the supported case; the caller turns
# it into the original positivity error.
struct VirtualGaugeUnsupported <: Exception
  message::String
end

# The averaged generator of one bosonic site on the Bargmann lattice |m) = a†^m|0>,
# e_{m,m'} = |m)(m'|:  L0 e_{m,m'} = λ(m,m') e_{m,m'} + loss m m' e_{m-1,m'-1}.
# `diagonal` lists the terms c (a†^k a^k) ρ (a†^k' a^k'), each of eigenvalue c m^(k) m'^(k')
# in falling factorials, and `loss` is the coefficient of a ρ a†.
struct BosonLadder
  diagonal::Vector{Tuple{Int,Int,LatticeScalar}}
  loss::LatticeScalar
end

falling_factorial(a::Int, l::Int) = a < l ? 0 : prod(a - t for t in 0:(l - 1); init=1)

function ladder_eigenvalue(ladder::BosonLadder, m::Int, m2::Int)
  value = zero(LatticeScalar)
  for (k, k2, c) in ladder.diagonal
    value += c * (falling_factorial(m, k) * falling_factorial(m2, k2))
  end
  return value
end

native_boson_ladder(::NativeRepresentation, _) = nothing

# `nothing` unless the representation is a single bosonic site and every term of L0 is either
# number-diagonal or the loss sandwich a ρ a†.
function native_boson_ladder(
  representation::AlgebraicLiouvilleRepresentation, L0::AlgebraSuperoperator
)
  sites = representation.algebra.sites
  (length(sites) == 1 && sites[1].kind == BOSON_SITE) || return nothing
  diagonal = Tuple{Int,Int,LatticeScalar}[]
  loss = zero(LatticeScalar)
  for ((X, M), c) in L0.terms
    if X == [0, 1] && M == [1, 0]
      loss += LatticeScalar(c)
    elseif X[1] == X[2] && M[1] == M[2]
      push!(diagonal, (X[1], M[1], LatticeScalar(c)))
    else
      return nothing
    end
  end
  return BosonLadder(diagonal, loss)
end

function checked_gap(ladder::BosonLadder, m, m2, u, u2)
  gap = ladder_eigenvalue(ladder, m, m2) - ladder_eigenvalue(ladder, u, u2)
  iszero(gap) &&
    throw(VirtualGaugeUnsupported("the averaged generator is resonant on the Fock lattice"))
  return gap
end

# Right eigenvector components r_j on e_{m-j,m'-j}, j = 0:min(m,m').
function lattice_right_eigenvector(ladder::BosonLadder, m::Int, m2::Int)
  r = LatticeScalar[one(LatticeScalar)]
  for j in 1:min(m, m2)
    push!(
      r,
      ladder.loss * ((m - j + 1) * (m2 - j + 1)) * r[end] /
      checked_gap(ladder, m, m2, m - j, m2 - j),
    )
  end
  return r
end

# Left eigenvector components l_i on e*_{m+i,m'+i}, i = 0:imax.
function lattice_left_eigenvector(ladder::BosonLadder, m::Int, m2::Int, imax::Int)
  l = LatticeScalar[one(LatticeScalar)]
  for i in 1:imax
    push!(
      l,
      ladder.loss * ((m + i) * (m2 + i)) * l[end] /
      checked_gap(ladder, m, m2, m + i, m2 + i),
    )
  end
  return l
end

# Action of a one-site superoperator on e_{a,b}, as a dictionary over the image lattice point.
function lattice_action(Y::AlgebraSuperoperator, a::Int, b::Int)
  image = Dict{Tuple{Int,Int},LatticeScalar}()
  for ((X, M), c) in Y.terms
    k, l = X
    k2, l2 = M
    f = falling_factorial(a, l) * falling_factorial(b, k2)
    f == 0 && continue
    key = (a - l + k, b - k2 + l2)
    (key[1] < 0 || key[2] < 0) && continue
    image[key] = get(image, key, zero(LatticeScalar)) + LatticeScalar(c) * f
  end
  return image
end

# A term c (a†^k a^l) ρ (a†^k' a^l') moves e_{a,b} by (k - l, l' - k'), so the neutral
# component reaches at most k left rungs.
function lattice_shift_bound(Y::AlgebraSuperoperator)
  bound = 0
  for ((X, M), _) in Y.terms
    bound = max(bound, X[1], M[2])
  end
  return bound
end

# The lattice invariants ℓ_{m,m'}(Y) = <L_{m,m'}|Y|R_{m,m'}> of several superoperators at one
# lattice point. They vanish on the range of ad L0 and are the resonance conditions of
# [L0, S] = Y on the Fock lattice.
function lattice_invariants(
  ladder::BosonLadder, Ys::AbstractVector{<:AlgebraSuperoperator}, m::Int, m2::Int
)
  bound = maximum(lattice_shift_bound, Ys; init=0)
  r = lattice_right_eigenvector(ladder, m, m2)
  l = lattice_left_eigenvector(ladder, m, m2, bound)
  values = zeros(LatticeScalar, length(Ys))
  for (index, Y) in pairs(Ys), (j, rj) in pairs(r)
    for ((a, b), c) in lattice_action(Y, m - (j - 1), m2 - (j - 1))
      i = a - m
      (i >= 0 && b - m2 == i && i <= bound) || continue
      values[index] += rj * l[i + 1] * c
    end
  end
  return values
end

function lattice_invariant(ladder::BosonLadder, Y::AlgebraSuperoperator, m::Int, m2::Int)
  return only(lattice_invariants(ladder, [Y], m, m2))
end

function max_lattice_invariant(ladder::BosonLadder, Y::AlgebraSuperoperator, window::Int)
  worst = zero(Rational{BigInt})
  for m in 0:window, m2 in 0:window
    v = lattice_invariant(ladder, Y, m, m2)
    worst = max(worst, abs(real(v)), abs(imag(v)))
  end
  return worst
end

struct VirtualParameter
  kind::Symbol
  first::Int
  second::Int
end

function virtual_unit(::Type{T}, n::Int, p::VirtualParameter) where {T}
  C = exact_zeros(T, n, n)
  if p.kind == :diagonal
    C[p.first, p.first] = one(T)
  elseif p.kind == :real
    C[p.first, p.second] = one(T)
    C[p.second, p.first] = one(T)
  elseif p.kind == :imag
    C[p.first, p.second] = im
    C[p.second, p.first] = -im
  end
  return C
end

function virtual_superoperator(
  representation::AlgebraicLiouvilleRepresentation{T,R}, p::VirtualParameter
) where {T,R}
  n = length(representation.frame)
  algebra = representation.algebra
  if p.kind == :hamiltonian
    H = algebra_operator(algebra, [([p.first, p.first], one(T))])
    return native_gksl(representation, H, exact_zeros(T, n, n))
  end
  H = AlgebraOperator{T,R}(algebra, Dict{Monomial,T}())
  return native_gksl(representation, H, virtual_unit(T, n, p))
end

# Neutral polynomial superoperators of degree at most `degree`: Hamiltonians a†^k a^k and
# Hermitian Kossakowski units between frame monomials of equal charge.
function neutral_virtual_parameters(
  representation::AlgebraicLiouvilleRepresentation, small::Vector{Int}, degree::Int
)
  frame = representation.frame
  parameters = VirtualParameter[
    VirtualParameter(:hamiltonian, k, 0) for k in 1:(degree ÷ 2)
  ]
  charge(μ) = frame[μ][1] - frame[μ][2]
  for (i, μ) in pairs(small), ν in small[i:end]
    charge(μ) == charge(ν) || continue
    if μ == ν
      push!(parameters, VirtualParameter(:diagonal, μ, μ))
    else
      push!(parameters, VirtualParameter(:real, μ, ν))
      push!(parameters, VirtualParameter(:imag, μ, ν))
    end
  end
  return parameters
end

function invariant_matrix(ladder::BosonLadder, Ys::Vector{<:AlgebraSuperoperator}, window)
  rows = Vector{Rational{BigInt}}[]
  for m in 0:window, m2 in 0:window
    values = lattice_invariants(ladder, Ys, m, m2)
    push!(rows, real.(values))
    push!(rows, imag.(values))
  end
  return Rational{BigInt}[row[j] for row in rows, j in eachindex(Ys)]
end

function virtual_gauge_combination(
  ::Type{T}, parameters::Vector{VirtualParameter}, x, n::Int
) where {T}
  C = exact_zeros(T, n, n)
  for (k, p) in pairs(parameters)
    (iszero(x[k]) || p.kind == :hamiltonian) && continue
    C += T(x[k]) * virtual_unit(T, n, p)
  end
  return C
end

# Terminal virtual gauge. At the last retained order the static gauge S_N is never stored. The
# defect Y = E_N - V̂_N is a polynomial GKSL superoperator whose lattice invariants vanish, so
# [L0, S_N] = Y has a (formal) solution on the Fock lattice and E_N has the same
# quasi-spectrum as V̂_N. Charged sectors lie in the range of ad L0 and are cancelled; the
# neutral polynomial parameters are restricted to the nullspace of the invariants.
function native_virtual_static_step(
  representation::AlgebraicLiouvilleRepresentation{T,R},
  ladder::BosonLadder,
  L0,
  residual,
  known::AbstractMatrix{<:Number},
  active::AbstractMatrix{<:Number},
  weights::AbstractVector{<:Number},
  order::Int,
) where {T,R}
  algebra = representation.algebra
  frame = representation.frame
  n = length(frame)
  charge(μ) = frame[μ][1] - frame[μ][2]
  small = Int[μ for μ in 1:n if sum(frame[μ]) <= order]

  parameters = neutral_virtual_parameters(representation, small, order)
  superoperators = AlgebraSuperoperator{T,R}[
    virtual_superoperator(representation, p) for p in parameters
  ]
  null = exact_nullspace(invariant_matrix(ladder, superoperators, 2 * order))
  images = Matrix{T}[virtual_gauge_combination(T, parameters, x, n) for x in null]

  # Charged Kossakowski entries are cancelled; the cross units between the charges the
  # active channels mix stay free, because the dark quotient ties them to neutral entries.
  charges = sort!(
    unique(charge(μ) for j in axes(active, 2) for μ in 1:n if !iszero(active[μ, j]))
  )
  CV = native_kossakowski(representation, residual)
  target = T[charge(μ) == charge(ν) ? CV[μ, ν] : known[μ, ν] for μ in 1:n, ν in 1:n]
  for μ in small, ν in small
    (charge(μ) in charges && charge(ν) in charges && charge(μ) < charge(ν)) || continue
    E = exact_zeros(T, n, n)
    E[μ, ν] = one(T)
    E[ν, μ] = one(T)
    push!(images, E)
    F = exact_zeros(T, n, n)
    F[μ, ν] = im
    F[ν, μ] = -im
    push!(images, F)
  end

  p = length(images)
  solution = native_exact_static_solve(
    target,
    known,
    active,
    R[real(w) for w in weights],
    images,
    representation.metric,
    Matrix{R}(LinearAlgebra.I, p, p),
  )

  coordinates = zeros(Rational{BigInt}, length(parameters))
  for (k, x) in pairs(null)
    coordinates .+= Rational{BigInt}(solution.coordinates[k]) .* x
  end
  H = AlgebraOperator{T,R}(algebra, Dict{Monomial,T}())
  for (k, q) in pairs(parameters)
    if q.kind == :hamiltonian && !iszero(coordinates[k])
      H += algebra_operator(algebra, [([q.first, q.first], T(coordinates[k]))])
    end
  end
  defect = native_gksl(representation, H, solution.coefficient - CV)
  max_lattice_invariant(ladder, defect, 5 * order) == 0 || throw(
    VirtualGaugeUnsupported(
      "the lattice invariants of the virtual defect do not vanish beyond the fitted window",
    ),
  )

  E = residual + defect
  Hamiltonian = native_hamiltonian(representation, E)
  native_kossakowski(representation, E) == solution.coefficient || throw(
    ArgumentError("virtual static step failed to reconstruct its Kossakowski coefficient")
  )
  E == native_gksl(representation, Hamiltonian, solution.coefficient) || throw(
    ArgumentError("virtual static step failed its Hamiltonian/Kossakowski reconstruction")
  )
  return NativeStaticStep(
    E,
    zero(L0),
    Hamiltonian,
    solution,
    solution.correction,
    T.(solution.newborn_weights),
    true,
    defect,
  )
end

function virtual_gauge_message(error::NativePositivityError, reason::AbstractString)
  return NativePositivityError(
    error.message * "; the terminal virtual gauge is not available: " * reason
  )
end

# Static step of one order. Polynomial gauge families are tried first. When all of them fail
# at the last retained order, the virtual gauge is tried; below the last order a virtual
# gauge cannot be used because later orders need S_n as an operator.
function native_order_step(
  representation::NativeRepresentation,
  inverse::HomologicalInverse,
  families::AbstractVector{<:NativeGaugeFamily},
  L0,
  residual,
  known::AbstractMatrix{<:Number},
  active::AbstractMatrix{<:Number},
  weights::AbstractVector{<:Number},
  tol::Real,
  order::Int,
  last_order::Int,
)
  try
    return native_static_step(
      representation, inverse, families, L0, residual, known, active, weights, tol
    )
  catch error
    error isa NativePositivityError || rethrow()
    ladder = native_boson_ladder(representation, L0)
    ladder === nothing && throw(
      virtual_gauge_message(
        error,
        "it needs a single bosonic site with number-diagonal averaged generator plus loss",
      ),
    )
    order == last_order || throw(
      ArgumentError(
        "order $order of the GKSL normal form needs a static gauge that is not polynomial, " *
        "and a virtual gauge is possible only at the last retained order $last_order. " *
        "An explicit static gauge would be required; lower the requested order.",
      ),
    )
    try
      return native_virtual_static_step(
        representation, ladder, L0, residual, known, active, weights, order
      )
    catch virtual_error
      virtual_error isa Union{VirtualGaugeUnsupported,NativePositivityError} || rethrow()
      throw(virtual_gauge_message(error, virtual_error.message))
    end
  end
end

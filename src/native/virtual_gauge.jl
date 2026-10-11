const VIRTUAL_FIT_WINDOW = 2
const VIRTUAL_VERIFICATION_WINDOW = 5
const LatticeScalar = Complex{Rational{BigInt}}

struct VirtualGaugeUnsupported <: Exception
  message::String
end

struct BosonLadder
  diagonal::Vector{Tuple{Int,Int,LatticeScalar}}
  loss::LatticeScalar
end

function boson_ladder(
  representation::AlgebraicLiouvilleRepresentation, L0::AlgebraSuperoperator
)
  sites = representation.algebra.sites
  length(sites) == 1 || throw(
    VirtualGaugeUnsupported(
      "it needs a single bosonic site, but the model has $(length(sites))"
    ),
  )
  sites[1].kind == BOSON_SITE || throw(VirtualGaugeUnsupported("it needs a bosonic site"))
  diagonal = Tuple{Int,Int,LatticeScalar}[]
  loss = zero(LatticeScalar)
  for ((X, M), c) in L0.terms
    if X == [0, 1] && M == [1, 0]
      loss += LatticeScalar(c)
    elseif X[1] == X[2] && M[1] == M[2]
      push!(diagonal, (X[1], M[1], LatticeScalar(c)))
    else
      throw(
        VirtualGaugeUnsupported(
          "it needs a number-diagonal averaged generator plus loss, " *
          "but the averaged generator has the term $X ρ $M",
        ),
      )
    end
  end
  return BosonLadder(diagonal, loss)
end

falling_factorial(a::Int, l::Int) = a < l ? 0 : prod(a - t for t in 0:(l - 1); init=1)

function ladder_eigenvalue(ladder::BosonLadder, m::Int, m2::Int)
  value = zero(LatticeScalar)
  for (k, k2, c) in ladder.diagonal
    value += c * (falling_factorial(m, k) * falling_factorial(m2, k2))
  end
  return value
end

function check_nonresonant(ladder::BosonLadder, window::Int)
  seen = Dict{LatticeScalar,Tuple{Int,Int}}()
  for m in 0:window, m2 in 0:window
    value = ladder_eigenvalue(ladder, m, m2)
    other = get!(seen, value, (m, m2))
    other == (m, m2) || throw(
      VirtualGaugeUnsupported(
        "the averaged generator is resonant on the Fock lattice: the points " *
        "$((m, m2)) and $other have the same eigenvalue",
      ),
    )
  end
  return true
end

function checked_gap(ladder::BosonLadder, point, other)
  gap = ladder_eigenvalue(ladder, point...) - ladder_eigenvalue(ladder, other...)
  iszero(gap) && throw(
    VirtualGaugeUnsupported(
      "the averaged generator is resonant on the Fock lattice: the points " *
      "$point and $other have the same eigenvalue",
    ),
  )
  return gap
end

function lattice_right_eigenvector(ladder::BosonLadder, m::Int, m2::Int)
  r = LatticeScalar[one(LatticeScalar)]
  for j in 1:min(m, m2)
    step = ladder.loss * ((m - j + 1) * (m2 - j + 1)) * r[end]
    push!(r, step / checked_gap(ladder, (m, m2), (m - j, m2 - j)))
  end
  return r
end

function lattice_left_eigenvector(ladder::BosonLadder, m::Int, m2::Int, bound::Int)
  l = LatticeScalar[one(LatticeScalar)]
  for i in 1:bound
    step = ladder.loss * ((m + i) * (m2 + i)) * l[end]
    push!(l, step / checked_gap(ladder, (m, m2), (m + i, m2 + i)))
  end
  return l
end

const NeutralTerm = Tuple{Int,Int,Int,LatticeScalar}

function neutral_terms(Y::AlgebraSuperoperator)
  terms = NeutralTerm[]
  for (((k, l), (k2, l2)), c) in Y.terms
    k - l == l2 - k2 && push!(terms, (k - l, l, k2, LatticeScalar(c)))
  end
  return terms
end

function lattice_invariants(
  ladder::BosonLadder, Ys::Vector{Vector{NeutralTerm}}, m::Int, m2::Int
)
  bound = maximum(Y -> maximum(first, Y; init=0), Ys; init=0)
  r = lattice_right_eigenvector(ladder, m, m2)
  l = lattice_left_eigenvector(ladder, m, m2, bound)
  values = zeros(LatticeScalar, length(Ys))
  for (index, Y) in pairs(Ys), (s, l1, k2, c) in Y, j in 0:min(s, m, m2)
    values[index] +=
      c *
      r[j + 1] *
      l[s - j + 1] *
      (falling_factorial(m - j, l1) * falling_factorial(m2 - j, k2))
  end
  return values
end

function invariant_matrix(ladder::BosonLadder, Ys::Vector{Vector{NeutralTerm}}, window::Int)
  rows = Vector{Rational{BigInt}}[]
  for m in 0:window, m2 in 0:window
    values = lattice_invariants(ladder, Ys, m, m2)
    push!(rows, real.(values))
    push!(rows, imag.(values))
  end
  return Rational{BigInt}[row[j] for row in rows, j in eachindex(Ys)]
end

function check_invariants(ladder::BosonLadder, Y::AlgebraSuperoperator, window::Int)
  terms = [neutral_terms(Y)]
  for m in 0:window, m2 in 0:window
    iszero(only(lattice_invariants(ladder, terms, m, m2))) || throw(
      VirtualGaugeUnsupported(
        "the lattice invariant of the defect at $((m, m2)) does not vanish"
      ),
    )
  end
  return true
end

function hermitian_units(::Type{T}, n::Int, μ::Int, ν::Int) where {T}
  real_unit = exact_zeros(T, n, n)
  real_unit[μ, ν] = one(T)
  real_unit[ν, μ] = one(T)
  μ == ν && return [real_unit]
  imag_unit = exact_zeros(T, n, n)
  imag_unit[μ, ν] = im
  imag_unit[ν, μ] = -im
  return [real_unit, imag_unit]
end

function native_virtual_static_step(
  representation::AlgebraicLiouvilleRepresentation{T,R},
  L0::AlgebraSuperoperator,
  residual,
  known::AbstractMatrix{<:Number},
  active::AbstractMatrix{<:Number},
  weights::AbstractVector{<:Number},
  order::Int,
) where {T,R}
  ladder = boson_ladder(representation, L0)
  check_nonresonant(ladder, VIRTUAL_VERIFICATION_WINDOW * order)
  algebra = representation.algebra
  frame = representation.frame
  n = length(frame)
  charge(μ) = frame[μ][1] - frame[μ][2]
  small = Int[μ for μ in 1:n if sum(frame[μ]) <= order]
  no_frame = exact_zeros(T, n, n)
  no_hamiltonian = AlgebraOperator{T,R}(algebra, Dict{Monomial,T}())

  candidates = Tuple{AlgebraOperator{T,R},Matrix{T}}[
    (algebra_operator(algebra, [([k, k], one(T))]), no_frame) for k in 1:(order ÷ 2)
  ]
  for (i, μ) in pairs(small), ν in small[i:end]
    charge(μ) == charge(ν) || continue
    append!(candidates, [(no_hamiltonian, C) for C in hermitian_units(T, n, μ, ν)])
  end
  neutral = [neutral_terms(native_gksl(representation, H, C)) for (H, C) in candidates]
  null = exact_nullspace(invariant_matrix(ladder, neutral, VIRTUAL_FIT_WINDOW * order))
  images = Matrix{T}[
    sum(T(x[k]) * candidates[k][2] for k in eachindex(x); init=no_frame) for x in null
  ]

  charges = sort!(
    unique(charge(μ) for j in axes(active, 2) for μ in 1:n if !iszero(active[μ, j]))
  )
  for μ in small, ν in small
    (charge(μ) < charge(ν) && charge(μ) in charges && charge(ν) in charges) || continue
    append!(images, hermitian_units(T, n, μ, ν))
  end

  CV = native_kossakowski(representation, residual)
  target = T[charge(μ) == charge(ν) ? CV[μ, ν] : known[μ, ν] for μ in 1:n, ν in 1:n]
  p = length(images)
  solution = native_exact_static_solve(
    target,
    Matrix{T}(known),
    Matrix{T}(active),
    R[real(w) for w in weights],
    images,
    representation.metric,
    Matrix{R}(LinearAlgebra.I, p, p),
  )

  H = no_hamiltonian
  for (coordinate, x) in zip(solution.coordinates, null), k in 1:(order ÷ 2)
    iszero(x[k]) && continue
    H += scale_operator(T(coordinate * x[k]), candidates[k][1])
  end
  defect = native_gksl(representation, H, solution.coefficient - CV)
  check_invariants(ladder, defect, VIRTUAL_VERIFICATION_WINDOW * order)

  E = residual + defect
  Hamiltonian = native_hamiltonian(representation, E)
  native_kossakowski(representation, E) == solution.coefficient || throw(
    ArgumentError("virtual static step failed to reconstruct its Kossakowski coefficient")
  )
  E == native_gksl(representation, Hamiltonian, solution.coefficient) || throw(
    ArgumentError("virtual static step failed its Hamiltonian/Kossakowski reconstruction")
  )
  return NativeStaticStep(E, zero(L0), Hamiltonian, solution, true)
end

function native_order_step(
  representation::NativeRepresentation,
  families::AbstractVector{<:NativeGaugeFamily},
  L0,
  residual,
  known::AbstractMatrix{<:Number},
  active::AbstractMatrix{<:Number},
  weights::AbstractVector{<:Number},
  order::Int,
  last_order::Int,
)
  try
    return native_static_step(
      representation, families, L0, residual, known, active, weights
    )
  catch error
    error isa NativePositivityError || rethrow()
    return native_terminal_step(
      representation, L0, residual, known, active, weights, order, last_order
    )
  end
end

function native_terminal_step(
  representation::NativeRepresentation,
  L0,
  residual,
  known,
  active,
  weights,
  order::Int,
  last_order::Int,
)
  try
    order == last_order || throw(
      VirtualGaugeUnsupported(
        "a virtual gauge exists only at the last retained order $last_order"
      ),
    )
    return native_virtual_static_step(
      representation, L0, residual, known, active, weights, order
    )
  catch error
    error isa Union{VirtualGaugeUnsupported,NativePositivityError} || rethrow()
    throw(
      ArgumentError(
        "retained order $order of the GKSL normal form has no static gauge of " *
        "polynomial degree at most $NATIVE_GAUGE_DEGREE, and the virtual gauge does " *
        "not apply: $(error.message); lower the requested order",
      ),
    )
  end
end

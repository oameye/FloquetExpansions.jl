"""
    FloquetExpansion

Result of [`floquet_expansion`](@ref). Read it with [`effective_hamiltonian`](@ref) and
[`kick_operator`](@ref) rather than by field access.
"""
struct FloquetExpansion{G<:Gauge}
  H::PeriodicOperator
  K::Vector{PeriodicOperator}
  Kdot::Vector{PeriodicOperator}
  dressedH::Vector{PeriodicOperator}
  dressedKdot::Vector{PeriodicOperator}
  Heff::Vector{SQA.QAdd}
  order::Int
end

_tri(n::Int, j::Int) = (n * (n + 1)) ÷ 2 + j + 1

_weightH(j::Int) = im^j * (1 // factorial(j))
_weightKdot(j::Int) = im^j * (1 // factorial(j + 1))

function _dressed_node(
  K::Vector{PeriodicOperator}, prev, n::Int, j::Int, H::PeriodicOperator
)
  acc = zero(H)
  for k in 1:(n - j + 1)
    acc = acc + SQA.commutator(K[k], prev(k, j))
  end
  return acc
end

function _assemble_resolvent(
  dressedH::Vector{PeriodicOperator},
  dressedKdot::Vector{PeriodicOperator},
  n::Int,
  H::PeriodicOperator,
)
  R = zero(H)
  for j in 0:n
    R = R + _weightH(j) * dressedH[_tri(n, j)]
  end
  for j in 1:n
    R = R - _weightKdot(j) * dressedKdot[_tri(n, j)]
  end
  return R
end

function Base.show(io::IO, ::MIME"text/plain", vv::FloquetExpansion{G}) where {G}
  return print(io, "FloquetExpansion{", nameof(G), "} of order ", vv.order)
end

_reattach(X::SQA.QAdd, wd, n::Int) = iszero(n) ? X : wd^(-n) * X
_reattach(X::PeriodicOperator, n::Int) = iszero(n) ? X : X.wd^(-n) * X

"""
    floquet_expansion(H::PeriodicOperator, gauge::Gauge, order::Int) -> FloquetExpansion
    floquet_expansion(H::QAdd, w, t, gauge::Gauge, order::Int) -> FloquetExpansion

Expand a periodically driven Hamiltonian into a time-independent effective Hamiltonian and a
periodic kick operator, to the given `order` in `1/H.wd`. The drive frequency is bound to `H`.

The second form parses the time dependence first, so `H` may be given as an ordinary
time-dependent operator in the drive frequency `w` and time `t`.

Truncation follows the spec: `order = 1` is the rotating-wave approximation, and the error is
``\\mathcal{O}(w_d^{-\\text{order}})``. Note that the Liouvillian literature counts this shifted
by one.

The series is ASYMPTOTIC, not convergent. Beyond a problem-dependent optimal order, a higher
`order` makes the answer worse rather than better.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables w::Real t::Real g::Real;

julia> H = harmonics(w * (a' * a) + g * cos(w * t) * (a + a'), w, t);

julia> vv = floquet_expansion(H, VanVleck(), 1)
FloquetExpansion{VanVleck} of order 1

julia> effective_hamiltonian(vv)
w * a' * a
```

See also [`effective_hamiltonian`](@ref), [`kick_operator`](@ref), [`harmonics`](@ref).
"""
function floquet_expansion(H::PeriodicOperator, gauge::Gauge, order::Int)
  order >= 1 || throw(ArgumentError("order must be >= 1, got $(order)"))
  LinearAlgebra.ishermitian(H) || throw(
    ArgumentError(
      "the drive is not Hermitian: it must satisfy H_{-m} = H_m' (eq:fourierH)"
    ),
  )

  nodes = (order * (order + 1)) ÷ 2
  dressedH = [zero(H) for _ in 1:nodes]
  dressedKdot = [zero(H) for _ in 1:nodes]
  K = PeriodicOperator[]
  Kdot = PeriodicOperator[]
  Heff = SQA.QAdd[]

  for n in 0:(order - 1)
    dressedH[_tri(n, 0)] = n == 0 ? H : zero(H)

    for j in 1:n
      dressedH[_tri(n, j)] = _dressed_node(
        K, (k, _) -> dressedH[_tri(n - k, j - 1)], n, j, H
      )
    end

    for j in 1:n
      dressedKdot[_tri(n, j)] = _dressed_node(
        K, (k, j_) -> j_ == 1 ? Kdot[n - k + 1] : dressedKdot[_tri(n - k, j_ - 1)], n, j, H
      )
    end

    R = _assemble_resolvent(dressedH, dressedKdot, n, H)
    R = SQA.simplify(R)

    Heffn = SQA.simplify(time_average(R))
    push!(Heff, Heffn)

    if n < order - 1
      Knext = SQA.simplify(antiderivative(R - Heffn, gauge))
      push!(K, Knext)
      push!(Kdot, derivative(Knext))
    end
  end

  return FloquetExpansion{typeof(gauge)}(H, K, Kdot, dressedH, dressedKdot, Heff, order)
end

function floquet_expansion(H::SQA.QAdd, w, t, gauge::Gauge, order::Int)
  return floquet_expansion(harmonics(H, w, t), gauge, order)
end

"""
    effective_hamiltonian(vv::FloquetExpansion) -> QAdd
    effective_hamiltonian(vv::FloquetExpansion, n::Int) -> QAdd

The time-independent effective Hamiltonian ``H_\\text{eff}^{[N]} = \\sum_{k<N} H_\\text{eff}^{(k)}``,
or with `n` the order-`n` contribution alone. The drive frequency stored in `vv.H.wd` is
reattached, so the order-`n` piece carries ``w_d^{-n}``.
"""
function effective_hamiltonian(vv::FloquetExpansion)
  out = zero(SQA.QAdd)
  for n in 0:(vv.order - 1)
    out = out + _reattach(vv.Heff[n + 1], vv.H.wd, n)
  end
  return SQA.simplify(out)
end

function effective_hamiltonian(vv::FloquetExpansion, n::Int)
  0 <= n < vv.order ||
    throw(ArgumentError("order $(n) is outside 0:$(vv.order - 1) for this expansion"))
  return SQA.simplify(_reattach(vv.Heff[n + 1], vv.H.wd, n))
end

"""
    kick_operator(vv::FloquetExpansion) -> PeriodicOperator
    kick_operator(vv::FloquetExpansion, t) -> QAdd

The periodic kick operator ``K^{[N]} = \\sum_{k<N} K^{(k)}`` generating the micromotion, as
harmonics or, given a time variable `t`, as the time-dependent operator ``K(t)``.

Under [`VanVleck`](@ref) this has vanishing time average, which is what makes
[`effective_hamiltonian`](@ref) independent of the initial phase of the drive.
"""
function kick_operator(vv::FloquetExpansion)
  out = zero(vv.H)
  for k in 1:length(vv.K)
    out = out + _reattach(vv.K[k], k)
  end
  return out
end

# Deliberately no `kick_operator(vv, n::Int)`: it would capture `kick_operator(vv, 0)`, which a
# caller means as t = 0, and silently return an order instead of a time.
kick_operator(vv::FloquetExpansion, t) = kick_operator(vv)(t)

# THROWAWAY PROTOTYPE — do not import from production code.
#
# This prototype answers one design question: can the existing homological
# expansion operate on one periodic-generator container for both Hamiltonian
# and symbolic Liouvillian components?

module LiouvillianAPIPrototype

import SecondQuantizedAlgebra as SQA
using Symbolics: Num, @variables

abstract type Gauge end

struct VanVleck <: Gauge end

qadd(x::SQA.QAdd) = x
qadd(x::SQA.QField) = 1 * x

const ActionKey = Tuple{SQA.QAdd, SQA.QAdd}
const TermDict = Dict{ActionKey, SQA.CNum}

"""A sum of elementary left/right actions on a density operator.

The key `(A, B)` denotes the map ``rho -> A*rho*B``. SQA owns all
Hilbert-space operator algebra; this wrapper only collects equal maps and
composes them.
"""
struct Liouvillian
  terms::TermDict
end

function add_term!(terms::TermDict, key::ActionKey, coefficient::SQA.CNum)
  iszero(coefficient) && return terms
  updated = get(terms, key, convert(SQA.CNum, 0)) + coefficient
  iszero(updated) ? delete!(terms, key) : (terms[key] = updated)
  return terms
end

function raw_liouvillian(terms::TermDict)
  normalized = TermDict()
  sizehint!(normalized, length(terms))
  for (key, coefficient) in terms
    add_term!(normalized, key, coefficient)
  end
  return Liouvillian(normalized)
end

function action(left::SQA.QField, right::SQA.QField, coefficient = 1)
  return raw_liouvillian(TermDict(
    (qadd(left), qadd(right)) => convert(SQA.CNum, coefficient),
  ))
end

"""Construct the symbolic action ``-i[H, \\cdot]``."""
function hamiltonian_action(H::SQA.QField)
  Hq = qadd(H)
  identity = one(Hq)
  return action(Hq, identity, -im) + action(identity, Hq, im)
end

"""Construct the symbolic dissipator ``D[L]``."""
function dissipator(L::SQA.QField)
  Lq = qadd(L)
  identity = one(Lq)
  norm = adjoint(Lq) * Lq
  return action(Lq, adjoint(Lq)) +
         action(norm, identity, -1 // 2) +
         action(identity, norm, -1 // 2)
end

"""Construct ``L = -i[H, \\cdot] + Σ_C D[C] + Σ_j γ_j D[J_j]``."""
function Liouvillian(H::SQA.QField; collapse_operators = (), jumps = (), rates = ())
  length(jumps) == length(rates) ||
    throw(ArgumentError("jumps and rates must have equal lengths"))

  generator = hamiltonian_action(H)
  for collapse in collapse_operators
    generator = generator + dissipator(collapse)
  end
  for (jump, rate) in zip(jumps, rates)
    generator = generator + rate * dissipator(jump)
  end
  return generator
end

Base.iszero(L::Liouvillian) = isempty(L.terms)
Base.isempty(L::Liouvillian) = isempty(L.terms)

Base.zero(::Liouvillian) = Liouvillian(TermDict())
Base.zero(::Type{Liouvillian}) = Liouvillian(TermDict())

Base.:(==)(L::Liouvillian, R::Liouvillian) = L.terms == R.terms
Base.isequal(L::Liouvillian, R::Liouvillian) = isequal(L.terms, R.terms)
Base.hash(L::Liouvillian, h::UInt) = hash(:Liouvillian, hash(L.terms, h))

function Base.:+(L::Liouvillian, R::Liouvillian)
  terms = copy(L.terms)
  for (key, coefficient) in R.terms
    add_term!(terms, key, coefficient)
  end
  return Liouvillian(terms)
end

Base.:-(L::Liouvillian) = -1 * L
Base.:-(L::Liouvillian, R::Liouvillian) = L + (-R)

function scale(coefficient, L::Liouvillian)
  return scale(convert(SQA.CNum, coefficient), L)
end

function scale(coefficient::SQA.CNum, L::Liouvillian)
  iszero(L) && return zero(L)
  iszero(coefficient) && return zero(L)
  return raw_liouvillian(Dict(
    key => coefficient * value for (key, value) in L.terms
  ))
end

Base.:*(coefficient, L::Liouvillian) = scale(coefficient, L)
Base.:*(L::Liouvillian, coefficient) = scale(coefficient, L)

"""Compose maps: `compose(A, B)(ρ) = A(B(ρ))`; `B` acts first."""
function compose(A::Liouvillian, B::Liouvillian)
  iszero(A) && return zero(A)
  iszero(B) && return zero(B)
  terms = TermDict()
  for ((left_A, right_A), coefficient_A) in A.terms,
      ((left_B, right_B), coefficient_B) in B.terms
    key = (left_A * left_B, right_B * right_A)
    add_term!(terms, key, coefficient_A * coefficient_B)
  end
  return Liouvillian(terms)
end

SQA.commutator(A::Liouvillian, B::Liouvillian) = compose(A, B) - compose(B, A)

function SQA.simplify(L::Liouvillian)
  terms = TermDict()
  for ((left, right), coefficient) in L.terms
    key = (SQA.simplify(left), SQA.simplify(right))
    add_term!(terms, key, coefficient)
  end
  return Liouvillian(terms)
end

Base.show(io::IO, L::Liouvillian) = print(io, isempty(L) ? "0" : "Liouvillian($(length(L.terms)) terms)")

"""One Fourier series type shared by Hamiltonian and Liouvillian generators."""
struct PeriodicGenerator{T}
  components::Dict{Int,T}
  wd::Num
  zero_component::T
end

function periodic_generator(components::AbstractDict{Int,T}, wd, zero_component::T) where {T}
  isconcretetype(T) || throw(ArgumentError("generator components must have a concrete element type"))
  kept = Dict{Int,T}()
  sizehint!(kept, length(components))
  for (harmonic, component) in components
    iszero(component) || (kept[harmonic] = component)
  end
  return PeriodicGenerator{T}(kept, Num(wd), zero_component)
end

function PeriodicGenerator(components::AbstractDict{Int,T}, wd) where {T}
  isempty(components) &&
    throw(ArgumentError("an empty PeriodicGenerator needs a zero component prototype"))
  return periodic_generator(components, wd, zero(first(values(components))))
end

function PeriodicGenerator(components::AbstractDict{Int,<:SQA.QField}, wd)
  normalized = Dict{Int,SQA.QAdd}(harmonic => qadd(component) for (harmonic, component) in components)
  isempty(normalized) &&
    throw(ArgumentError("an empty PeriodicGenerator needs a zero component prototype"))
  return periodic_generator(normalized, wd, zero(first(values(normalized))))
end

Base.getindex(G::PeriodicGenerator, harmonic::Int) = get(G.components, harmonic, G.zero_component)
Base.keys(G::PeriodicGenerator) = keys(G.components)
Base.length(G::PeriodicGenerator) = length(G.components)
Base.iszero(G::PeriodicGenerator) = isempty(G.components)
Base.isempty(G::PeriodicGenerator) = isempty(G.components)

function Base.:(==)(L::PeriodicGenerator, R::PeriodicGenerator)
  return isequal(L.wd, R.wd) && L.components == R.components && L.zero_component == R.zero_component
end

function Base.isequal(L::PeriodicGenerator, R::PeriodicGenerator)
  return isequal(L.wd, R.wd) && isequal(L.components, R.components) && isequal(L.zero_component, R.zero_component)
end

function Base.hash(G::PeriodicGenerator, h::UInt)
  return hash(:PeriodicGenerator, hash(G.wd, hash(G.components, hash(G.zero_component, h))))
end

function support(G::PeriodicGenerator)
  isempty(G) && return 0:-1
  harmonics = keys(G.components)
  return minimum(harmonics):maximum(harmonics)
end

function Base.show(io::IO, ::MIME"text/plain", G::PeriodicGenerator)
  if isempty(G)
    print(io, "PeriodicGenerator (zero)")
    return nothing
  end
  harmonics = sort!(collect(keys(G.components)))
  print(io, "PeriodicGenerator with harmonics ", first(harmonics), ":", last(harmonics))
  for harmonic in harmonics
    print(io, "\n  l = ", harmonic, "  =>  ", G.components[harmonic])
  end
  return nothing
end

Base.show(io::IO, G::PeriodicGenerator) = show(io, MIME"text/plain"(), G)

function Base.zero(G::PeriodicGenerator{T}) where {T}
  return periodic_generator(Dict{Int,T}(), G.wd, G.zero_component)
end

function check_frequency(L::PeriodicGenerator, R::PeriodicGenerator)
  isequal(L.wd, R.wd) || throw(ArgumentError("periodic generators use different drive frequencies"))
  return nothing
end

function Base.:+(L::PeriodicGenerator{T}, R::PeriodicGenerator{T}) where {T}
  check_frequency(L, R)
  harmonics = Dict{Int,T}(L.components)
  for (harmonic, component) in R.components
    harmonics[harmonic] = haskey(harmonics, harmonic) ? harmonics[harmonic] + component : component
  end
  return periodic_generator(harmonics, L.wd, L.zero_component)
end

Base.:-(L::PeriodicGenerator) = -1 * L
Base.:-(L::PeriodicGenerator{T}, R::PeriodicGenerator{T}) where {T} = L + (-R)

function Base.:*(coefficient, L::PeriodicGenerator{T}) where {T}
  return periodic_generator(
    Dict{Int,T}(harmonic => coefficient * component for (harmonic, component) in L.components),
    L.wd,
    coefficient * L.zero_component,
  )
end
Base.:*(L::PeriodicGenerator, coefficient) = coefficient * L

function SQA.commutator(L::PeriodicGenerator{T}, R::PeriodicGenerator{T}) where {T}
  check_frequency(L, R)
  out = Dict{Int,T}()
  for (left_harmonic, left) in L.components, (right_harmonic, right) in R.components
    harmonic = left_harmonic + right_harmonic
    term = SQA.commutator(left, right)
    out[harmonic] = haskey(out, harmonic) ? out[harmonic] + term : term
  end
  return periodic_generator(out, L.wd, L.zero_component)
end

time_average(G::PeriodicGenerator) = G[0]

function remove_average(G::PeriodicGenerator{T}) where {T}
  return periodic_generator(
    Dict{Int,T}(harmonic => component for (harmonic, component) in G.components if harmonic != 0),
    G.wd,
    G.zero_component,
  )
end

function derivative(G::PeriodicGenerator{T}) where {T}
  return periodic_generator(
    Dict{Int,T}(harmonic => (-im * harmonic) * component for (harmonic, component) in G.components),
    G.wd,
    G.zero_component,
  )
end

function antiderivative(G::PeriodicGenerator{T}, ::VanVleck) where {T}
  haskey(G.components, 0) && throw(ArgumentError("antiderivative requires zero average"))
  return periodic_generator(
    Dict{Int,T}(harmonic => (im // harmonic) * component for (harmonic, component) in G.components),
    G.wd,
    G.zero_component,
  )
end

function SQA.simplify(G::PeriodicGenerator{T}) where {T}
  return periodic_generator(
    Dict{Int,T}(harmonic => SQA.simplify(component) for (harmonic, component) in G.components),
    G.wd,
    G.zero_component,
  )
end

struct FloquetExpansion{G,P,E}
  generator::P
  kick_components::Vector{P}
  kick_derivative_components::Vector{P}
  dressed_generator::Vector{P}
  dressed_kick_derivative::Vector{P}
  effective_components::Vector{E}
  gauge::G
  order::Int
end

triindex(n::Int, j::Int) = (n * (n + 1)) ÷ 2 + j + 1
weight_generator(j::Int) = im^j * (1 // factorial(j))
weight_kick_derivative(j::Int) = im^j * (1 // factorial(j + 1))

function dressed_node(K::Vector{P}, previous, n::Int, j::Int, generator::P) where {P}
  result = zero(generator)
  for k in 1:(n - j + 1)
    result = result + SQA.commutator(K[k], previous(k, j))
  end
  return result
end

function assemble_resolvent(
  dressed_generator::Vector{P},
  dressed_kick_derivative::Vector{P},
  n::Int,
  generator::P,
) where {P}
  result = zero(generator)
  for j in 0:n
    result = result + weight_generator(j) * dressed_generator[triindex(n, j)]
  end
  for j in 1:n
    result = result - weight_kick_derivative(j) * dressed_kick_derivative[triindex(n, j)]
  end
  return result
end

function floquet_expansion(generator::P, gauge::G, order::Int) where {P<:PeriodicGenerator,G<:Gauge}
  order >= 1 || throw(ArgumentError("order must be >= 1"))

  nodes = (order * (order + 1)) ÷ 2
  dressed_generator = [zero(generator) for _ in 1:nodes]
  dressed_kick_derivative = [zero(generator) for _ in 1:nodes]
  K = P[]
  Kdot = P[]
  E = typeof(time_average(generator))
  effective = E[]

  for n in 0:(order - 1)
    dressed_generator[triindex(n, 0)] = n == 0 ? generator : zero(generator)

    for j in 1:n
      dressed_generator[triindex(n, j)] = dressed_node(
        K,
        (k, _) -> dressed_generator[triindex(n - k, j - 1)],
        n,
        j,
        generator,
      )
    end

    for j in 1:n
      dressed_kick_derivative[triindex(n, j)] = dressed_node(
        K,
        (k, j_) -> j_ == 1 ? Kdot[n - k + 1] : dressed_kick_derivative[triindex(n - k, j_ - 1)],
        n,
        j,
        generator,
      )
    end

    resolvent = SQA.simplify(assemble_resolvent(dressed_generator, dressed_kick_derivative, n, generator))
    effective_n = SQA.simplify(time_average(resolvent))
    push!(effective, effective_n)

    if n < order - 1
      next_kick = SQA.simplify(antiderivative(remove_average(resolvent), gauge))
      push!(K, next_kick)
      push!(Kdot, derivative(next_kick))
    end
  end

  return FloquetExpansion(
    generator,
    K,
    Kdot,
    dressed_generator,
    dressed_kick_derivative,
    effective,
    gauge,
    order,
  )
end

reattach(component, wd, n::Int) = iszero(n) ? component : wd^(-n) * component

function effective_generator(expansion::FloquetExpansion)
  result = zero(expansion.effective_components[1])
  for n in 0:(expansion.order - 1)
    result = result + reattach(expansion.effective_components[n + 1], expansion.generator.wd, n)
  end
  return SQA.simplify(result)
end

function effective_generator(expansion::FloquetExpansion, n::Int)
  0 <= n < expansion.order ||
    throw(ArgumentError("order $(n) is outside 0:$(expansion.order - 1)"))
  return SQA.simplify(reattach(expansion.effective_components[n + 1], expansion.generator.wd, n))
end

function micromotion(expansion::FloquetExpansion)
  result = zero(expansion.generator)
  for (order, kick) in enumerate(expansion.kick_components)
    result = result + reattach(kick, expansion.generator.wd, order)
  end
  return result
end

function micromotion(expansion::FloquetExpansion, n::Int)
  1 <= n < expansion.order ||
    throw(ArgumentError("order $(n) is outside 1:$(expansion.order - 1)"))
  return SQA.simplify(reattach(expansion.kick_components[n], expansion.generator.wd, n))
end

function demo()
  space = SQA.FockSpace(:prototype)
  a = SQA.Destroy(space, :a)
  @variables ω::Real γ::Real

  L0 = Liouvillian(a' * a; jumps = (a,), rates = (γ,))
  L1 = Liouvillian(a)
  generator = PeriodicGenerator(Dict(0 => L0, 1 => L1, -1 => L1), ω)
  expansion = floquet_expansion(generator, VanVleck(), 2)

  println("Liouvillian type: ", typeof(L0))
  println("Periodic type: ", typeof(generator))
  println("Commutator: ", SQA.commutator(L0, L1))
  println("Effective type: ", typeof(effective_generator(expansion)))
  println("Effective generator: ", effective_generator(expansion))
  println("Micromotion type: ", typeof(micromotion(expansion)))
  return expansion
end

end

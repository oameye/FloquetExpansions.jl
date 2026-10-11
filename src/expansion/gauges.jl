"""
    Gauge

Supertype of gauges that fix the free integration constant in [`antiderivative`](@ref).
"""
abstract type Gauge end

"""
    ExpansionAlgorithm

Abstract selector for the algorithm that computes a Floquet expansion within a [`Gauge`](@ref).

The gauge fixes which effective generator and micromotion [`floquet_expansion`](@ref) returns;
the algorithm fixes only how they are computed, so algorithms for the same gauge return the same
retained coefficients. They do differ in cost. Pass one through the gauge constructor, as in
`VanVleck(; algorithm=BlochFeshbach())`. The implemented algorithms are [`HoriDeprit`](@ref) and
[`BlochFeshbach`](@ref).
"""
abstract type ExpansionAlgorithm end

"""
    HoriDeprit <: ExpansionAlgorithm

Select the Hori–Deprit Lie-transform algorithm, which solves for the micromotion generator order by
order through nested commutators of the drive harmonics.

`HoriDeprit()` is the default algorithm of [`VanVleck`](@ref), so `VanVleck()` and
`VanVleck(; algorithm=HoriDeprit())` select the same gauge. The construction is described in
[High-frequency expansion](@ref high-frequency-expansion-theory).

See also [`BlochFeshbach`](@ref).
"""
struct HoriDeprit <: ExpansionAlgorithm end

"""
    BlochFeshbach <: ExpansionAlgorithm

Select the Bloch/Feshbach wave-operator algorithm. It solves Bloch's equation for the wave operator
order by order and normalizes the result to the van Vleck micromotion and effective generator.

For the same generator and `order`, the retained [`effective_component`](@ref) and
[`micromotion`](@ref) coefficients equal those of [`HoriDeprit`](@ref). The construction is
described in [Bloch/Feshbach projection](@ref bloch-feshbach-theory).

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables ω::Real t::Real Δ::Real g::Real;

julia> H = harmonics(Δ * (a' * a) + g * cos(ω * t) * (a' * a' + a * a), ω, t);

julia> bf = floquet_expansion(H, VanVleck(; algorithm=BlochFeshbach()), 3);

julia> hd = floquet_expansion(H, VanVleck(), 3);

julia> iszero(simplify(effective_generator(bf) - effective_generator(hd)))
true
```
"""
struct BlochFeshbach <: ExpansionAlgorithm end

"""
    VanVleck(; algorithm=HoriDeprit())

Select the van Vleck gauge, ``\\langle \\mathcal{K} \\rangle = 0``. This gives the micromotion
generator zero period average and makes the effective generator independent of the drive's
initial phase. The condition fixes every retained effective and micromotion coefficient uniquely.

`algorithm` selects the [`ExpansionAlgorithm`](@ref) that computes these coefficients without
changing them: [`HoriDeprit`](@ref), the default, or [`BlochFeshbach`](@ref).

# References

The van Vleck construction follows [VanVleck1929](@cite), and its Floquet-space formulation
follows [Rahav2003](@cite), [Eckardt2015](@cite), and [Bukov2015](@cite).
"""
struct VanVleck{A<:ExpansionAlgorithm} <: Gauge
  algorithm::A
end

VanVleck(; algorithm::ExpansionAlgorithm=HoriDeprit()) = VanVleck(algorithm)

"""
    GKSLNormalForm(; algorithm=HoriDeprit(), negative=(), verification=8)

Select the GKSL normal-form gauge for Liouvillian expansions. It adds a static similarity to the
van Vleck micromotion, chosen order by order so that every truncation of the effective generator
is of GKSL form and therefore completely positive. The effective generator keeps the van Vleck
spectrum through the retained order, but its retained coefficients differ from those of
[`VanVleck`](@ref) wherever the static gauge is nonzero. Without dissipation the two coincide.

`algorithm` selects the [`ExpansionAlgorithm`](@ref) that computes these coefficients without
changing them: [`HoriDeprit`](@ref), the default, or [`BlochFeshbach`](@ref).

For a Liouvillian the result is read with [`effective_generator`](@ref), [`hamiltonian`](@ref),
[`channels`](@ref) and [`kossakowski`](@ref), as after [`positive_completion`](@ref). For a
Hamiltonian it holds the van Vleck coefficients and no completed model, so `channels` and
`kossakowski` are not defined on it. The construction is described in
[GKSL normal form](@ref gksl-normal-form-theory).

# Symbolic parameters

Every free symbol in the coefficients is a parameter, and the result is symbolic in them. It is
reconstructed rather than derived: the exact numeric expansion is evaluated at rational points
and every coefficient is rebuilt as an exact rational function of the parameters. The result is
claimed on one declared region only, where every parameter is positive except those listed in
`negative`, which are negative, and the drive frequency is positive.

- **Checked.** After the fit, every coefficient is compared exactly with the numeric expansion at
  `verification` fresh rational points of the region. Each point scales every parameter
  independently by a factor between ``1/27`` and ``27`` of a value between ``1/9`` and
  ``9/2``. At each point the frame, the channel onsets, and the orders carrying a static gauge
  must equal those of the reference point, or an `ArgumentError` is thrown.
  `factorization(expansion)` returns `(; parameters, signs, points)`, the verified points
  in the order of `parameters`, so the claim can be re-checked by substituting a point
  into the result and into the numeric expansion.
- **Not guaranteed.** Agreement at any other point. A solver branch confined to a part of the
  region that no point samples goes unseen. The sampled points therefore bound the claim.
- **Reported.** [`positivity_conditions`](@ref) lists the signed parameters, the drive
  frequency, and every channel rate including its power of the drive frequency. The rates
  are nonnegative only under these conditions, which is why the drive frequency appears.

Every coefficient must be a polynomial of moderate degree in the parameters, with integer
powers; a coefficient that is rational in the parameters throws an `ArgumentError`. At most six
parameters are supported and the cost grows quickly with their number, so substitute exact
values for the parameters that need not stay free.

# Restrictions

- Supported subsystems are bosonic modes, position and momentum, `NLevelSpace` and `PauliSpace`
  levels, and spins of any size. Bosonic modes are not truncated.
- Operators on one Hilbert space must share one name; `Position(h, :x)` and `Momentum(h, :p)`
  together throw an `ArgumentError`.
- A numeric coefficient must be an exact `Integer` or `Rational`. A floating-point coefficient
  throws an `ArgumentError`, so write `4 // 5` rather than `0.8`.
- `order` must be at least `1`.
- An order whose static gauge is not polynomial throws an `ArgumentError`; lower the requested
  order.
- [`positive_completion`](@ref) throws on the result, which is already of GKSL form.

# References

The truncation problem that this gauge removes is discussed in [Schnell2021](@cite) and
[Ikeda2021](@cite).

# Examples

A driven Kerr oscillator with loss. The drive induces a second channel at order four:

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables ω::Real t::Real;

julia> H = (1 // 2) * a' * a + (3 // 10) * a' * a * a' * a + (2 // 5) * (a + a') * cos(ω * t);

julia> x = floquet_expansion(H, ω, t, GKSLNormalForm(), 5; channels=(jump(a, 4 // 5),));

julia> length(channels(x))
2

julia> channels(x)[2]
((72//3125) / (ω^4))𝒟[a * a - 2 * a' * a]
```

A driven qubit with decay and free parameters. The drive induces a channel whose rate is
reconstructed symbolically, and the verified points are returned by [`factorization`](@ref):

```jldoctest
julia> q = PauliSpace(:qubit); σx = Pauli(q, :σ, 1); σy = Pauli(q, :σ, 2); σz = Pauli(q, :σ, 3);

julia> @variables ω::Real t::Real Δ::Real F::Real κ::Real;

julia> H = Δ * σz / 2 + F * cos(ω * t) * σx;

julia> x = floquet_expansion(H, ω, t, GKSLNormalForm(), 4; channels=(jump((σx - im * σy) / 2, κ),));

julia> channels(x)[2]
((2(F^2)*κ) / (ω^2))𝒟[1//2 + 1//2 * σz]

julia> length(factorization(x).points)
8
```
"""
struct GKSLNormalForm{A<:ExpansionAlgorithm} <: Gauge
  algorithm::A
  negative::Vector{Symbolics.Num}
  verification::Int
end

function GKSLNormalForm(;
  algorithm::ExpansionAlgorithm=HoriDeprit(), negative=Symbolics.Num[], verification::Int=8
)
  verification >= 1 || throw(ArgumentError("verification must be >= 1"))
  return GKSLNormalForm(algorithm, Symbolics.Num[negative...], verification)
end

"""
    antiderivative(X::PeriodicGenerator, gauge::Gauge) -> PeriodicGenerator

Integrate `X` with respect to dimensionless drive time, with `gauge` fixing the free
integration constant:

```math
(\\partial_t^{-1} X)_l = \\frac{i}{l} X_l \\quad (l \\neq 0)
```

`X` must have vanishing time average; pass `X - time_average(X)` if that is not already
true.

# Notes

Weights remain exact rationals, so an `OverflowError` from `Rational{Int}` is possible at
high order rather than a silent loss of precision.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables ω::Real t::Real;

julia> X = harmonics(a * expim(2ω * t), ω, t);

julia> derivative(antiderivative(X, VanVleck())) == X
true
```
"""
function antiderivative(G::PeriodicGenerator{T}, ::VanVleck) where {T}
  haskey(G.components, 0) && throw(ArgumentError("antiderivative requires zero average"))
  return periodic_generator(
    Dict{Int,T}(
      harmonic => (im // harmonic) * component for (harmonic, component) in G.components
    ),
    G.wd,
    G.zero_component,
  )
end

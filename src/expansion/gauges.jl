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
    GKSLNormalForm(; algorithm=HoriDeprit(), negative=())

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

Coefficients may be polynomials in at most eight symbolic parameters with exact rational
coefficients. The expansion is then computed exactly over the field of rational functions in
those parameters, and the channel rates are rational functions. Every parameter and the drive
frequency are taken positive, except those listed in `negative`, which are taken negative. The
[`positivity_conditions`](@ref) of the result list these signs and the channel rates.

Every sign the construction decides must follow from these signs alone, so the result holds
throughout the sign region. A sign that does not, such as that of `κ - γ`, throws an
`ArgumentError`. The result further assumes generic parameters: no expression the expansion
divides by vanishes. At an exceptional point substitute the values before expanding.

# Restrictions

- Supported subsystems are bosonic modes, position and momentum, `NLevelSpace` and `PauliSpace`
  levels, and spins of any size. Bosonic modes are not truncated.
- Operators on one Hilbert space must share one name; `Position(h, :x)` and `Momentum(h, :p)`
  together throw an `ArgumentError`.
- Every coefficient must be exact: a floating-point coefficient throws an `ArgumentError`, so
  write `4 // 5` rather than `0.8`. A symbolic coefficient must be a polynomial with integer
  powers of its parameters.
- `order` must be at least `1`.
- An order with no static gauge of polynomial degree at most three throws an `ArgumentError`;
  lower the requested order.
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

The same model with symbolic parameters gives the rate of that channel exactly:

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables ω::Real t::Real Δ::Real χ::Real F::Real κ::Real;

julia> H = Δ * a' * a + χ * a' * a * a' * a + F * (a + a') * cos(ω * t);

julia> x = floquet_expansion(H, ω, t, GKSLNormalForm(), 5; channels=(jump(a, κ),));

julia> channels(x)[2]
((2(F^2)*κ*(χ^2)) / (ω^4))𝒟[a * a - 2 * a' * a]
```
"""
struct GKSLNormalForm{A<:ExpansionAlgorithm} <: Gauge
  algorithm::A
  negative::Vector{Symbolics.Num}
end

function GKSLNormalForm(; algorithm::ExpansionAlgorithm=HoriDeprit(), negative=())
  return GKSLNormalForm(algorithm, Symbolics.Num[negative...])
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

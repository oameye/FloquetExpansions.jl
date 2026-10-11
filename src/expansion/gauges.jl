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
    GKSLNormalForm(; algorithm=HoriDeprit())

Select the GKSL normal-form gauge for Liouvillian expansions. The effective Liouvillian is of
GKSL form, and therefore completely positive, at every order. The oscillatory part of the
micromotion is fixed as in [`VanVleck`](@ref); its static part is a similarity
``W = e^{S}``, ``S = \\sum_n S_n/ω^n``, chosen order by order so that every truncation of the
effective generator is of GKSL form.

The effective generator is ``W^{-1} \\mathcal{L}_{\\mathrm{VV}} W``. It has the van Vleck spectrum
through the retained order and coincides with the van Vleck generator wherever every ``S_n``
vanishes, which is always the case without dissipation. The finite generator is assembled from
graded jump amplitudes, so it equals the retained series up to a positive remainder beyond the
retained order. The result is read with [`effective_generator`](@ref), [`hamiltonian`](@ref),
[`channels`](@ref) and [`kossakowski`](@ref).

The expansion runs in an exact operator algebra without truncation: bosonic modes, position
and momentum, `NLevelSpace` and `PauliSpace` levels, and spins as the size-independent
``\\mathfrak{su}(2)`` algebra. `algorithm` selects [`HoriDeprit`](@ref), the default, or
[`BlochFeshbach`](@ref); both return the same coefficients.

# Restrictions

- Every coefficient must be an exact `Integer` or `Rational` and numeric. A floating-point
  coefficient or a symbolic parameter other than the drive frequency throws an `ArgumentError`,
  so write `4 // 5` rather than `0.8`.
- `order` must be at least `1`.
- An order whose static gauge is not polynomial throws an `ArgumentError`; lower the requested
  order.
- [`positive_completion`](@ref) throws on the result, which is already of GKSL form.

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
"""
struct GKSLNormalForm{A<:ExpansionAlgorithm} <: Gauge
  algorithm::A
end

GKSLNormalForm(; algorithm::ExpansionAlgorithm=HoriDeprit()) = GKSLNormalForm(algorithm)

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

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

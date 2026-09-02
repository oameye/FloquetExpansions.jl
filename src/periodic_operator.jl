"""
    Gauge

Supertype of the gauge choices that fix the free constant in [`antiderivative`](@ref).
"""
abstract type Gauge end

"""
    VanVleck()

The van Vleck gauge, ``\\langle K \\rangle = 0``: the micromotion generator has vanishing period
average, which makes the effective generator independent of the initial phase of the drive.
"""
struct VanVleck <: Gauge end

function validate_component_type(::Type{T}) where {T}
  isconcretetype(T) ||
    throw(ArgumentError("generator components must have a concrete element type"))
  return nothing
end

"""
    PeriodicGenerator(components::AbstractDict{Int, T}, wd)
    PeriodicGenerator(components::AbstractDict{Int, T}, wd, zero_component::T) where {T}

A ``T``-periodic generator held by Fourier harmonic, with drive frequency `wd`,

```math
G(t) = \\sum_l G_l \\, e^{-i l \\omega_d t}.
```

Index `l` selects a harmonic. A missing harmonic and a zero harmonic are the same thing: both
are absent from the stored components, so `G[l]` returns zero for any `l` and `keys(G)` lists
only the nonzero ones. The frequency is part of `G`, so generators can only be combined when they
use the same Fourier basis.

For an empty generator, the optional `zero_component` prototype fixes the component type. The
component type determines the algebra used by addition, commutators, differentiation, and
simplification.

`G` is Hermitian exactly when `G[-l] == G[l]'` for every `l`, which is what `ishermitian` checks
for Hamiltonian components and what `adjoint` produces.

# Examples

```jldoctest
julia> using LinearAlgebra: ishermitian

julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables w::Real t::Real;

julia> G = harmonics(a * expim(w * t) + a' * expim(-w * t), w, t)
PeriodicGenerator with harmonics -1:1
  l = -1  =>  a
  l = 1  =>  a'

julia> G[1]
a'

julia> ishermitian(G)
true
```

See also [`time_average`](@ref), [`derivative`](@ref), [`antiderivative`](@ref).
"""
struct PeriodicGenerator{T}
  components::Dict{Int,T}
  wd::Symbolics.Num
  zero_component::T

  function PeriodicGenerator{T}(
    components::Dict{Int,T}, wd::Symbolics.Num, zero_component::T
  ) where {T}
    validate_component_type(T)
    kept = Dict{Int,T}()
    sizehint!(kept, length(components))
    for (harmonic, component) in components
      iszero(component) || (kept[harmonic] = component)
    end
    return new{T}(kept, wd, zero_component)
  end
end

@inline component_pairs(G::PeriodicGenerator) = pairs(G.components)

const PeriodicScalar = Union{Real,Complex,Symbolics.Num,SQA.CNum}

function periodic_generator(
  components::AbstractDict{Int,T}, wd::Symbolics.Num, zero_component::T
) where {T}
  return PeriodicGenerator{T}(Dict{Int,T}(components), wd, zero_component)
end

function PeriodicGenerator(
  components::AbstractDict{Int,T}, wd::Symbolics.Num, zero_component::T
) where {T}
  return periodic_generator(components, wd, zero_component)
end

function PeriodicGenerator(components::AbstractDict{Int,T}, wd::Symbolics.Num) where {T}
  validate_component_type(T)
  zero_component = isempty(components) ? zero(T) : zero(first(values(components)))
  return PeriodicGenerator(components, wd, zero_component)
end

function PeriodicGenerator(components::AbstractDict{Int,<:SQA.QField}, wd::Symbolics.Num)
  normalized = Dict{Int,SQA.QAdd}(
    harmonic => qadd(component) for (harmonic, component) in components
  )
  return periodic_generator(normalized, wd, zero(SQA.QAdd))
end

function Base.getindex(G::PeriodicGenerator, harmonic::Int)
  return get(G.components, harmonic, G.zero_component)
end
Base.keys(G::PeriodicGenerator) = keys(G.components)
Base.length(G::PeriodicGenerator) = length(G.components)
Base.iszero(G::PeriodicGenerator) = isempty(G.components)
Base.isempty(G::PeriodicGenerator) = isempty(G.components)

function Base.:(==)(L::PeriodicGenerator, R::PeriodicGenerator)
  return isequal(L.wd, R.wd) &&
         L.components == R.components &&
         L.zero_component == R.zero_component
end

function Base.isequal(L::PeriodicGenerator, R::PeriodicGenerator)
  return isequal(L.wd, R.wd) &&
         isequal(L.components, R.components) &&
         isequal(L.zero_component, R.zero_component)
end

function Base.hash(G::PeriodicGenerator, h::UInt)
  return hash(:PeriodicGenerator, hash(G.wd, hash(G.components, hash(G.zero_component, h))))
end

"""
    support(G::PeriodicGenerator) -> UnitRange{Int}

Smallest range of harmonic indices containing every nonzero harmonic of `G`, or `0:-1` if `G`
is zero.
"""
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

function check_frequency(L::PeriodicGenerator, R::PeriodicGenerator)
  isequal(L.wd, R.wd) ||
    throw(ArgumentError("periodic generators use different drive frequencies"))
  return nothing
end

function Base.:+(L::PeriodicGenerator{T}, R::PeriodicGenerator{T}) where {T}
  check_frequency(L, R)
  harmonics = Dict{Int,T}(L.components)
  for (harmonic, component) in R.components
    harmonics[harmonic] =
      haskey(harmonics, harmonic) ? harmonics[harmonic] + component : component
  end
  return periodic_generator(harmonics, L.wd, L.zero_component)
end

Base.:-(L::PeriodicGenerator) = -1 * L
Base.:-(L::PeriodicGenerator{T}, R::PeriodicGenerator{T}) where {T} = L + (-R)

function Base.:+(L::PeriodicGenerator{T}, component::T) where {T}
  return L + PeriodicGenerator(Dict(0 => component), L.wd)
end
Base.:+(component::T, L::PeriodicGenerator{T}) where {T} = L + component
Base.:-(L::PeriodicGenerator{T}, component::T) where {T} = L + (-component)
function Base.:-(component::T, L::PeriodicGenerator{T}) where {T}
  return PeriodicGenerator(Dict(0 => component), L.wd) - L
end

function Base.zero(G::PeriodicGenerator{T}) where {T}
  return periodic_generator(Dict{Int,T}(), G.wd, G.zero_component)
end

function Base.:*(coefficient::PeriodicScalar, L::PeriodicGenerator{T}) where {T}
  return periodic_generator(
    Dict{Int,T}(
      harmonic => coefficient * component for (harmonic, component) in L.components
    ),
    L.wd,
    coefficient * L.zero_component,
  )
end
Base.:*(L::PeriodicGenerator, coefficient::PeriodicScalar) = coefficient * L

"""
    commutator(K::PeriodicGenerator, X::PeriodicGenerator) -> PeriodicGenerator

Commutator of two periodic generators, which is the harmonic convolution

```math
[K, X]_l = \\sum_p [K_p,\\, X_{l-p}]
```

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables w::Real t::Real;

julia> K = harmonics(a * expim(-w * t), w, t); X = harmonics(a' * expim(-w * t), w, t);

julia> commutator(K, X)[2]
1
```
"""
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

"""
    time_average(G::PeriodicGenerator) -> T

Average of `G` over one period, which is its zeroth harmonic `G[0]`.
"""
time_average(G::PeriodicGenerator{T}) where {T} = G[0]

function remove_average(G::PeriodicGenerator{T}) where {T}
  return periodic_generator(
    Dict{Int,T}(
      harmonic => component for (harmonic, component) in G.components if harmonic != 0
    ),
    G.wd,
    G.zero_component,
  )
end

"""
    derivative(G::PeriodicGenerator) -> PeriodicGenerator

Time derivative in units of the drive frequency, `dG/dt` divided by ``\\omega_d``:

```math
(\\partial_t G)_l = -i l \\, G_l.
```
"""
function derivative(G::PeriodicGenerator{T}) where {T}
  return periodic_generator(
    Dict{Int,T}(
      harmonic => (-im * harmonic) * component for (harmonic, component) in G.components
    ),
    G.wd,
    G.zero_component,
  )
end

"""
    antiderivative(X::PeriodicGenerator, gauge::Gauge) -> PeriodicGenerator

Inverse of [`derivative`](@ref), with `gauge` fixing the free integration constant:

```math
(\\partial_t^{-1} X)_l = \\frac{i}{l} X_l \\quad (l \\neq 0)
```

`X` must have vanishing time average; pass `X - time_average(X)` if that is not already true.
Weights stay exact rationals, so an `OverflowError` from `Rational{Int}` is possible at high
order rather than a silent loss of precision.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables w::Real t::Real;

julia> X = harmonics(a * expim(2w * t), w, t);

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

"""
    simplify(G::PeriodicGenerator) -> PeriodicGenerator

Simplify every component of `G` symbolically.
"""
function SQA.simplify(G::PeriodicGenerator{T}) where {T}
  return periodic_generator(
    Dict{Int,T}(
      harmonic => SQA.simplify(component) for (harmonic, component) in G.components
    ),
    G.wd,
    G.zero_component,
  )
end

"""
    adjoint(G::PeriodicGenerator) -> PeriodicGenerator

Hermitian adjoint, harmonic by harmonic: `adjoint(G)[l] == adjoint(G[-l])`.
"""
function Base.adjoint(G::PeriodicGenerator{SQA.QAdd})
  return PeriodicGenerator(
    Dict{Int,SQA.QAdd}(
      -harmonic => adjoint(component) for (harmonic, component) in G.components
    ),
    G.wd,
    adjoint(G.zero_component),
  )
end

"""
    ishermitian(G::PeriodicGenerator) -> Bool

True when `G[-l] == adjoint(G[l])` for every harmonic, i.e. when the time-dependent generator
is Hermitian at every time for Hamiltonian components.
"""
function LinearAlgebra.ishermitian(G::PeriodicGenerator{SQA.QAdd})
  for (harmonic, component) in G.components
    iszero(SQA.simplify(adjoint(component) - G[-harmonic])) || return false
  end
  return true
end

"""
    (G::PeriodicGenerator)(t)

Rebuild the time-dependent generator from its Fourier components. The return
value has the component type of `G`.
"""
function (G::PeriodicGenerator)(t::Symbolics.Num)
  return sum(
    SQA.expim(-harmonic * G.wd * t) * component for (harmonic, component) in G.components;
    init=zero(G.zero_component),
  )
end

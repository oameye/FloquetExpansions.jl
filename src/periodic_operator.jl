"""
    Gauge

Supertype of the gauge choices that fix the free constant in [`antiderivative`](@ref).
"""
abstract type Gauge end

"""
    VanVleck()

The van Vleck gauge, ``\\langle K \\rangle = 0``: the kick operator has vanishing period
average, which makes the effective Hamiltonian independent of the initial phase of the drive.
"""
struct VanVleck <: Gauge end

"""
    PeriodicOperator(components::AbstractDict{Int, <:QField}, wd)

A ``T``-periodic operator held by Fourier harmonic, with drive frequency `wd`,

```math
X(t) = \\sum_l X_l \\, e^{-i l \\omega_d t}
```

Index `l` selects a harmonic. A missing harmonic and a zero harmonic are the same thing: both
are absent from the stored components, so `X[l]` returns zero for any `l` and `keys(X)` lists
only the nonzero ones. The frequency is part of `X`, so operators can only be combined when they
use the same Fourier basis.

`X` is Hermitian exactly when `X[-l] == X[l]'` for every `l`, which is what
`ishermitian` checks and what `adjoint` produces.

# Examples

```jldoctest
julia> using LinearAlgebra: ishermitian

julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables w::Real t::Real;

julia> X = harmonics(a * expim(w * t) + a' * expim(-w * t), w, t)
PeriodicOperator with harmonics -1:1
  l = -1  =>  a
  l =  1  =>  a'

julia> X[1]
a'

julia> ishermitian(X)
true
```

See also [`time_average`](@ref), [`derivative`](@ref), [`antiderivative`](@ref).
"""
struct PeriodicOperator
  components::Dict{Int,SQA.QAdd}
  wd::Symbolics.Num

  function PeriodicOperator(components::Dict{Int,SQA.QAdd}, wd)
    kept = Dict{Int,SQA.QAdd}()
    sizehint!(kept, length(components))
    for (l, Xl) in components
      iszero(Xl) || (kept[l] = Xl)
    end
    return new(kept, Symbolics.Num(wd))
  end
end

# `Op` is not `<: QAdd`, so every ingest path has to promote.
promote_qadd(x::SQA.QAdd) = x
promote_qadd(x::SQA.QSym) = 1 * x

function PeriodicOperator(components::AbstractDict{Int,<:SQA.QField}, wd)
  return PeriodicOperator(
    Dict{Int,SQA.QAdd}(l => promote_qadd(Xl) for (l, Xl) in components), wd
  )
end

Base.getindex(X::PeriodicOperator, l::Int) = get(X.components, l, zero(SQA.QAdd))
Base.keys(X::PeriodicOperator) = keys(X.components)
Base.length(X::PeriodicOperator) = length(X.components)
Base.iszero(X::PeriodicOperator) = isempty(X.components)
Base.isempty(X::PeriodicOperator) = isempty(X.components)

function Base.:(==)(X::PeriodicOperator, Y::PeriodicOperator)
  return isequal(X.wd, Y.wd) && X.components == Y.components
end
function Base.isequal(X::PeriodicOperator, Y::PeriodicOperator)
  return isequal(X.wd, Y.wd) && isequal(X.components, Y.components)
end
function Base.hash(X::PeriodicOperator, h::UInt)
  return hash(:PeriodicOperator, hash(X.wd, hash(X.components, h)))
end

"""
    support(X::PeriodicOperator) -> UnitRange{Int}

Smallest range of harmonic indices containing every nonzero harmonic of `X`, or `0:-1` if
`X` is zero.
"""
function support(X::PeriodicOperator)
  isempty(X.components) && return 0:-1
  return minimum(keys(X.components)):maximum(keys(X.components))
end

function Base.show(io::IO, ::MIME"text/plain", X::PeriodicOperator)
  if isempty(X.components)
    print(io, "PeriodicOperator (zero)")
    return nothing
  end
  ls = sort!(collect(keys(X.components)))
  print(io, "PeriodicOperator with harmonics ", first(ls), ":", last(ls))
  pad = maximum(length ∘ string, ls)
  for l in ls
    print(io, "\n  l = ", lpad(l, pad), "  =>  ", X.components[l])
  end
  return nothing
end

Base.show(io::IO, X::PeriodicOperator) = show(io, MIME"text/plain"(), X)

function check_frequency(X::PeriodicOperator, Y::PeriodicOperator)
  isequal(X.wd, Y.wd) || throw(
    ArgumentError(
      "cannot combine PeriodicOperators with different drive frequencies: " *
      "$(X.wd) and $(Y.wd)",
    ),
  )
  return nothing
end

function addto!(out::Dict{Int,SQA.QAdd}, l::Int, Xl::SQA.QAdd)
  out[l] = haskey(out, l) ? out[l] + Xl : Xl
  return out
end

function Base.:+(X::PeriodicOperator, Y::PeriodicOperator)
  check_frequency(X, Y)
  out = copy(X.components)
  for (l, Yl) in Y.components
    addto!(out, l, Yl)
  end
  return PeriodicOperator(out, X.wd)
end

function Base.:-(X::PeriodicOperator)
  return PeriodicOperator(Dict{Int,SQA.QAdd}(l => -Xl for (l, Xl) in X.components), X.wd)
end
Base.:-(X::PeriodicOperator, Y::PeriodicOperator) = X + (-Y)

# The static operand is a DC harmonic; `R - time_average(R)` is the recursion's own idiom.
function Base.:+(X::PeriodicOperator, c::SQA.QField)
  return X + PeriodicOperator(Dict{Int,SQA.QAdd}(0 => promote_qadd(c)), X.wd)
end
function Base.:+(c::SQA.QField, X::PeriodicOperator)
  return X + PeriodicOperator(Dict{Int,SQA.QAdd}(0 => promote_qadd(c)), X.wd)
end
function Base.:-(X::PeriodicOperator, c::SQA.QField)
  return X + PeriodicOperator(Dict{Int,SQA.QAdd}(0 => -promote_qadd(c)), X.wd)
end
function Base.:-(c::SQA.QField, X::PeriodicOperator)
  return PeriodicOperator(Dict{Int,SQA.QAdd}(0 => promote_qadd(c)), X.wd) + (-X)
end

function Base.:*(c::Number, X::PeriodicOperator)
  return PeriodicOperator(Dict{Int,SQA.QAdd}(l => c * Xl for (l, Xl) in X.components), X.wd)
end
Base.:*(X::PeriodicOperator, c::Number) = c * X

Base.zero(X::PeriodicOperator) = PeriodicOperator(Dict{Int,SQA.QAdd}(), X.wd)

"""
    adjoint(X::PeriodicOperator) -> PeriodicOperator

Hermitian adjoint, harmonic by harmonic: `X'[l] == X[-l]'`.
"""
function Base.adjoint(X::PeriodicOperator)
  return PeriodicOperator(
    Dict{Int,SQA.QAdd}(-l => adjoint(Xl) for (l, Xl) in X.components), X.wd
  )
end

"""
    ishermitian(X::PeriodicOperator) -> Bool

True when `X[-l] == X[l]'` for every harmonic, i.e. when `X(t)` is Hermitian at every `t`.
"""
function LinearAlgebra.ishermitian(X::PeriodicOperator)
  for (l, Xl) in X.components
    iszero(SQA.simplify(adjoint(Xl) - X[-l])) || return false
  end
  return true
end

"""
    simplify(X::PeriodicOperator) -> PeriodicOperator

Simplify every coefficient of `X` symbolically.
"""
function SQA.simplify(X::PeriodicOperator)
  return PeriodicOperator(
    Dict{Int,SQA.QAdd}(l => SQA.simplify(Xl) for (l, Xl) in X.components), X.wd
  )
end

"""
    commutator(K::PeriodicOperator, X::PeriodicOperator) -> PeriodicOperator

Commutator of two periodic operators, which is the harmonic convolution

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
function SQA.commutator(K::PeriodicOperator, X::PeriodicOperator)
  check_frequency(K, X)
  out = Dict{Int,SQA.QAdd}()
  for (p, Kp) in K.components, (q, Xq) in X.components
    c = SQA.commutator(Kp, Xq)
    iszero(c) || addto!(out, p + q, c)
  end
  return PeriodicOperator(out, K.wd)
end

"""
    derivative(X::PeriodicOperator) -> PeriodicOperator

Time derivative in units of the drive frequency, `dX/dt` divided by ``\\omega_d``:

```math
(\\partial_t X)_l = -i l \\, X_l
```
"""
function derivative(X::PeriodicOperator)
  return PeriodicOperator(
    Dict{Int,SQA.QAdd}(l => (-im * l) * Xl for (l, Xl) in X.components), X.wd
  )
end

"""
    time_average(X::PeriodicOperator) -> QAdd

Average of `X` over one period, which is its zeroth harmonic `X[0]`.
"""
time_average(X::PeriodicOperator) = X[0]

"""
    antiderivative(X::PeriodicOperator, gauge::Gauge) -> PeriodicOperator

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
function antiderivative(X::PeriodicOperator, ::VanVleck)
  haskey(X.components, 0) && throw(
    ArgumentError(
      "antiderivative requires a vanishing time average, but harmonic 0 is present; " *
      "pass `X - time_average(X)`",
    ),
  )
  return PeriodicOperator(
    Dict{Int,SQA.QAdd}(l => ((1 // l) * im) * Xl for (l, Xl) in X.components), X.wd
  )
end

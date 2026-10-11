"""
    CarrierEmbedding(a, frequencies; t)

Embed one physical bosonic mode `a` into carrier modes for quantum harmonic balance. Every
carrier frequency `ωⱼ` contributes a co-rotating mode `aⱼ₊` and a counter-rotating mode
`aⱼ₋`, and the physical mode is represented by

```math
a \\mapsto \\frac{1}{\\sqrt{λ_H}} \\sum_j \\left(e^{-iω_j t} a_{j+} + e^{+iω_j t} a_{j-}\\right),
\\qquad λ_H = 2N_h,
```

with ``N_h`` the number of carriers. The weight ``1/\\sqrt{λ_H}`` is kept exact, so the map
preserves ``[a, a^\\dagger] = 1`` and every product of the physical algebra.

Frequencies are symbolic expressions in shared symbols with exact coefficients. A relation
between carriers is written into the expressions themselves, as in
`[ω₁, (ω₁ + ω₂) / 2, ω₂]` for a midpoint carrier; [`harmonic_balance`](@ref) keeps exactly
the processes whose frequency then vanishes identically. Floating-point coefficients of a
symbol and frequencies that depend on `t` throw an `ArgumentError`.

The carrier modes live on a new product space and are returned by [`carrier_modes`](@ref).
The embedding acts on expressions in `a` alone; [`reconstruct`](@ref) maps a physical
operator into the carrier algebra.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables ω₁::Real ω₂::Real t::Real;

julia> embedding = CarrierEmbedding(a, [ω₁, ω₂]; t);

julia> length(carrier_modes(embedding))
2

julia> I = reconstruct(embedding, a);

julia> commutator(I, I')
1
```

See also [`harmonic_balance`](@ref), [`reconstruct`](@ref), [`carrier_modes`](@ref).
"""
struct CarrierEmbedding
  physical::SQA.Op
  frequencies::Vector{Symbolics.Num}
  modes::Vector{Tuple{SQA.Op,SQA.Op}}
  time::Symbolics.Num
  image::SQA.QAdd
end

function CarrierEmbedding(a::SQA.Op, frequencies::AbstractVector; t::Symbolics.Num)
  SQA.is_destroy(a) || throw(
    ArgumentError("a carrier embedding maps a bosonic annihilation operator, got $(a)")
  )
  isempty(frequencies) &&
    throw(ArgumentError("a carrier embedding needs at least one carrier frequency"))
  carriers = Symbolics.Num[validated_carrier_frequency(ω, t) for ω in frequencies]
  name = SQA.operator_name(a)
  count = length(carriers)
  spaces = [
    SQA.FockSpace(Symbol(name, j, sector)) for j in 1:count for sector in ("₊", "₋")
  ]
  space = reduce(SQA.:⊗, spaces)
  modes = Tuple{SQA.Op,SQA.Op}[
    (
      SQA.Destroy(space, Symbol(name, j, "₊"), 2j - 1),
      SQA.Destroy(space, Symbol(name, j, "₋"), 2j),
    ) for j in 1:count
  ]
  weight = sqrt(Symbolics.Num(1 // (2count)))
  image = sum(
    weight * (SQA.expim(-ω * t) * plus + SQA.expim(ω * t) * minus) for
    (ω, (plus, minus)) in zip(carriers, modes)
  )
  return CarrierEmbedding(a, carriers, modes, t, image)
end

function validated_carrier_frequency(ω, t::Symbolics.Num)::Symbolics.Num
  frequency = Symbolics.Num(ω)
  any(isequal(Symbolics.value(t)), Symbolics.get_variables(frequency)) &&
    throw(ArgumentError("carrier frequency $(frequency) depends on the time $(t)"))
  validate_exact_frequency(frequency)
  return frequency
end

function contains_inexact_number(x)::Bool
  value = Symbolics.unwrap_const(Symbolics.unwrap(x))
  value isa AbstractFloat && return !isinteger(value)
  value isa Complex &&
    return contains_inexact_number(real(value)) || contains_inexact_number(imag(value))
  value isa Number && return false
  Symbolics.iscall(value) || return false
  return any(contains_inexact_number, Symbolics.arguments(value))
end

function validate_exact_frequency(frequency::Symbolics.Num)
  expanded = Symbolics.expand(frequency)
  (!isempty(Symbolics.get_variables(expanded)) && contains_inexact_number(expanded)) &&
    throw(
      ArgumentError(
        "frequency $(frequency) has a floating-point coefficient; write it with exact " *
        "coefficients, such as (ω₁ + ω₂) / 2, so that resonances cancel exactly",
      ),
    )
  return nothing
end

"""
    carrier_modes(embedding::CarrierEmbedding) -> Vector{Tuple{Op,Op}}

Return the carrier modes of `embedding` as `(aⱼ₊, aⱼ₋)` pairs, in the order of its carrier
frequencies. `aⱼ₊` carries the phase ``e^{-iω_j t}`` and `aⱼ₋` the phase ``e^{+iω_j t}`` in
the image of the physical mode.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables ω::Real t::Real;

julia> only(carrier_modes(CarrierEmbedding(a, [ω]; t)))
(a1₊, a1₋)
```

See also [`CarrierEmbedding`](@ref).
"""
carrier_modes(embedding::CarrierEmbedding) = copy(embedding.modes)

"""
    reconstruct(embedding::CarrierEmbedding, A) -> QAdd

Map a physical operator `A` into the carrier algebra, ``A \\mapsto I_t(A)``. Expectation
values of physical observables are read from a carrier state ``ρ`` as
``⟨A⟩(t) = \\mathrm{Tr}[I_t(A)\\,ρ]``. The map is an algebra homomorphism, so
`reconstruct(embedding, A * B) == reconstruct(embedding, A) * reconstruct(embedding, B)`.

`A` may contain only the embedded mode and its adjoint; any other operator throws an
`ArgumentError`.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables ω::Real t::Real;

julia> embedding = CarrierEmbedding(a, [ω]; t);

julia> reconstruct(embedding, a' * a) == reconstruct(embedding, a') * reconstruct(embedding, a)
true
```

See also [`CarrierEmbedding`](@ref), [`harmonic_balance`](@ref).
"""
function reconstruct(embedding::CarrierEmbedding, A::SQA.QField)
  operator = qadd(A)
  physical = embedding.physical
  for op in SQA.get_operators(operator)
    (isequal(op, physical) || isequal(op, adjoint(physical))) || throw(
      ArgumentError(
        "operator $(op) is not the embedded mode $(physical); a carrier embedding acts " *
        "on expressions in one physical mode",
      ),
    )
  end
  return SQA.substitute(operator, Dict(physical => embedding.image))
end

reconstruct(embedding::CarrierEmbedding, A::Number) = qadd(A * one(SQA.QAdd))

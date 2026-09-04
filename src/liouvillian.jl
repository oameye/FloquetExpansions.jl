const LiouvillianAction = Tuple{SQA.QAdd,SQA.QAdd}
const LiouvillianTerms = Dict{LiouvillianAction,SQA.CNum}
const LiouvillianScalar = Union{Number,Symbolics.Num,SQA.CNum}
const LiouvillianChannelCollection = Union{Tuple,AbstractVector}

abstract type LiouvillianChannel end

struct CollapseChannel{O<:SQA.QField} <: LiouvillianChannel
  operator::O
end

struct RateWeightedJump{O<:SQA.QField} <: LiouvillianChannel
  operator::O
  rate::SQA.CNum
  assumption::NonnegativeRateAssumption
end

"""
    Liouvillian

A symbolic linear map on density operators represented as a collected sum of elementary
actions ``ρ ↦ AρB``. The operator factors are SQA expressions; scalar coefficients are
symbolic SQA coefficients.

Use [`liouvillian`](@ref), [`hamiltonian_action`](@ref), or [`dissipator`](@ref) to construct
Liouvillians.
"""
struct Liouvillian
  terms::LiouvillianTerms
end

@inline function add_term!(
  L::Liouvillian, left::SQA.QAdd, right::SQA.QAdd, coefficient::SQA.CNum
)
  (iszero(left) || iszero(right) || iszero(coefficient)) && return L
  key = (left, right)
  updated = get(L.terms, key, convert(SQA.CNum, 0)) + coefficient
  iszero(updated) ? delete!(L.terms, key) : (L.terms[key] = updated)
  return L
end

function raw_liouvillian(terms::LiouvillianTerms)
  normalized = Liouvillian(LiouvillianTerms())
  sizehint!(normalized.terms, length(terms))
  for ((left, right), coefficient) in terms
    add_term!(normalized, left, right, coefficient)
  end
  return normalized
end

@inline term_pairs(L::Liouvillian) = pairs(L.terms)

"""
    terms(L::Liouvillian)

Return an iterator over the elementary terms in `L`. Each item is a
`(left_operator, right_operator, coefficient)` tuple representing
`ρ ↦ coefficient * left_operator * ρ * right_operator`.

The iteration order is unspecified. The returned iterator is independent of the sparse
storage used by `L`.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> L = liouvillian(a' * a; channels=(collapse(a),));

julia> all(length(item) == 3 for item in terms(L))
true
```

See also [`compose`](@ref).
"""
function terms(L::Liouvillian)
  return ((left, right, coefficient) for ((left, right), coefficient) in term_pairs(L))
end

function action(left::SQA.QField, right::SQA.QField, coefficient::LiouvillianScalar=1)
  return raw_liouvillian(
    LiouvillianTerms((qadd(left), qadd(right)) => convert(SQA.CNum, coefficient))
  )
end

"""
    hamiltonian_action(H::QField) -> Liouvillian

Construct the coherent density-operator action
``ρ ↦ -i[H, ρ]`` for the Hamiltonian `H`.

See also [`dissipator`](@ref), [`Liouvillian`](@ref).
"""
function hamiltonian_action(H::SQA.QField)
  Hq = qadd(H)
  identity = one(Hq)
  return action(Hq, identity, -im) + action(identity, Hq, im)
end

"""
    dissipator(L::QField) -> Liouvillian

Construct the symbolic dissipator
``D[L](ρ) = LρL† - (L†Lρ + ρL†L)/2``.

`L` is a complete collapse or jump operator. Use [`jump`](@ref) when a separate scalar rate
is part of the channel.

See also [`collapse`](@ref), [`jump`](@ref).
"""
function dissipator(L::SQA.QField)
  Lq = qadd(L)
  identity = one(Lq)
  norm = adjoint(Lq) * Lq
  return action(Lq, adjoint(Lq)) +
         action(norm, identity, -1 // 2) +
         action(identity, norm, -1 // 2)
end

function _microscopic_provenance(channels::LiouvillianChannelCollection)
  collapse_operators = SQA.QAdd[]
  jump_operators = SQA.QAdd[]
  jump_rates = SQA.CNum[]
  rate_assumptions = NonnegativeRateAssumption[]
  order = DissipativeSeedRef[]

  for channel in channels
    if channel isa CollapseChannel
      push!(collapse_operators, qadd(channel.operator))
      push!(order, DissipativeSeedRef(COLLAPSE_SEED, length(collapse_operators)))
    elseif channel isa RateWeightedJump
      push!(jump_operators, qadd(channel.operator))
      push!(jump_rates, channel.rate)
      push!(rate_assumptions, channel.assumption)
      push!(order, DissipativeSeedRef(JUMP_SEED, length(jump_operators)))
    else
      throw(
        ArgumentError("channels must contain only `collapse(...)` and `jump(...)` values")
      )
    end
  end

  return MicroscopicProvenance(
    collapse_operators, jump_operators, jump_rates, rate_assumptions, order
  )
end

function _liouvillian_from_provenance(H::SQA.QField, provenance::MicroscopicProvenance)
  generator = hamiltonian_action(H)
  for operator in provenance.collapse_operators
    generator = generator + dissipator(operator)
  end
  for i in eachindex(provenance.jump_operators)
    generator =
      generator + provenance.jump_rates[i] * dissipator(provenance.jump_operators[i])
  end
  return generator
end

"""
    liouvillian(H::QField; channels=()) -> Liouvillian

Construct
``ρ ↦ -i[H, ρ] + Σₐ D[Cₐ](ρ) + Σᵦ γᵦD[Jᵦ](ρ)``.

`channels` is a tuple or vector of [`collapse`](@ref) and [`jump`](@ref) values. A
[`jump`](@ref) rate is a physical rate: it must be real and is assumed nonnegative when
symbolic.

# Arguments

- `H`: Symbolic Hamiltonian expression.
- `channels`: Tuple or vector of channel values to add to the coherent action.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables γ::Real;

julia> L = liouvillian(a' * a; channels=(jump(a, γ),));

julia> L == hamiltonian_action(a' * a) + γ * dissipator(a)
true
```
"""
function liouvillian(H::SQA.QField; channels::LiouvillianChannelCollection=())
  return _liouvillian_from_provenance(H, _microscopic_provenance(channels))
end

@inline _channel_liouvillian(channel::CollapseChannel) = dissipator(channel.operator)
@inline _channel_liouvillian(channel::RateWeightedJump) =
  channel.rate * dissipator(channel.operator)

"""
    collapse(operator::QField) -> CollapseChannel

Create a channel from a complete collapse operator `L`. When passed to a
[`liouvillian`](@ref), it contributes the [`dissipator`](@ref) of `operator` with unit weight:
``D[L](ρ) = LρL† - (L†Lρ + ρL†L)/2``.
Here, `L` denotes `operator`; any amplitude or phase belonging to the collapse process is included
in it.

See also [`jump`](@ref), [`liouvillian`](@ref).
"""
function collapse(operator::SQA.QField)
  return CollapseChannel(operator)
end

function _validated_jump_rate(rate::LiouvillianScalar)
  coefficient = convert(SQA.CNum, rate)
  value = SQA.to_num(coefficient)
  imaginary_part = Symbolics.simplify(imag(value))
  iszero(imaginary_part) ||
    throw(ArgumentError("jump rate must be provably real; got `$rate`"))

  real_part = Symbolics.simplify(real(value))
  unwrapped = Symbolics.value(real_part)
  if unwrapped isa Real && unwrapped < 0
    throw(ArgumentError("jump rate must be nonnegative; got `$rate`"))
  end
  return coefficient
end

"""
    jump(operator::QField, rate) -> RateWeightedJump

Create a physical rate-weighted channel from a bare jump operator `J` and a separate rate
``γ``. When passed to a [`liouvillian`](@ref), it contributes ``γ D[J](ρ)``.

The rate must be provably real. A provably negative numeric rate is rejected immediately.
A symbolic real rate or real symbolic expression is accepted as an explicit physical
assumption that the whole expression satisfies ``γ ≥ 0``; FloquetExpansions does not split a
symbolic expression into separately signed factors. Use [`collapse`](@ref) when the complete
channel amplitude is already folded into the operator.

See also [`collapse`](@ref), [`liouvillian`](@ref).
"""
function jump(operator::SQA.QField, rate::LiouvillianScalar)
  coefficient = _validated_jump_rate(rate)
  return RateWeightedJump(operator, coefficient, NonnegativeRateAssumption(coefficient))
end

Base.iszero(L::Liouvillian) = isempty(L.terms)
Base.isempty(L::Liouvillian) = isempty(L.terms)

Base.zero(::Liouvillian) = Liouvillian(LiouvillianTerms())
Base.zero(::Type{Liouvillian}) = Liouvillian(LiouvillianTerms())

Base.one(::Type{Liouvillian}) = action(one(SQA.QAdd), one(SQA.QAdd))
Base.one(::Liouvillian) = one(Liouvillian)

Base.:(==)(L::Liouvillian, R::Liouvillian) = L.terms == R.terms
Base.isequal(L::Liouvillian, R::Liouvillian) = isequal(L.terms, R.terms)
Base.hash(L::Liouvillian, h::UInt) = hash(:Liouvillian, hash(L.terms, h))

function Base.:+(L::Liouvillian, R::Liouvillian)
  result = Liouvillian(copy(L.terms))
  for ((left, right), coefficient) in term_pairs(R)
    add_term!(result, left, right, coefficient)
  end
  return result
end

Base.:-(L::Liouvillian) = -1 * L
Base.:-(L::Liouvillian, R::Liouvillian) = L + (-R)

function scale(coefficient::LiouvillianScalar, L::Liouvillian)
  return scale(convert(SQA.CNum, coefficient), L)
end

function scale(coefficient::SQA.CNum, L::Liouvillian)
  iszero(L) && return zero(L)
  iszero(coefficient) && return zero(L)
  return raw_liouvillian(
    LiouvillianTerms(key => coefficient * value for (key, value) in L.terms)
  )
end

Base.:*(coefficient::LiouvillianScalar, L::Liouvillian) = scale(coefficient, L)
Base.:*(L::Liouvillian, coefficient::LiouvillianScalar) = scale(coefficient, L)

"""
    compose(A::Liouvillian, B::Liouvillian) -> Liouvillian

Compose maps so that `B` acts first and `A` acts second. For elementary actions
``cₐ Aₗρ Aᵣ`` and ``cᵦ Bₗρ Bᵣ``, the composed action is
``cₐcᵦ AₗBₗρ BᵣAᵣ``.

See also [`terms`](@ref).
"""
function compose(A::Liouvillian, B::Liouvillian)
  iszero(A) && return zero(A)
  iszero(B) && return zero(B)
  result = Liouvillian(LiouvillianTerms())
  for ((left_A, right_A), coefficient_A) in term_pairs(A),
    ((left_B, right_B), coefficient_B) in term_pairs(B)

    add_term!(result, left_A * left_B, right_B * right_A, coefficient_A * coefficient_B)
  end
  return result
end

SQA.commutator(A::Liouvillian, B::Liouvillian) = compose(A, B) - compose(B, A)

function SQA.simplify(L::Liouvillian)
  result = Liouvillian(LiouvillianTerms())
  for ((left, right), coefficient) in term_pairs(L)
    add_term!(result, SQA.simplify(left), SQA.simplify(right), coefficient)
  end
  return result
end

"""
    harmonics(L::Liouvillian, ω, t) -> PeriodicGenerator{Liouvillian}

Decompose a symbolic time-dependent Liouvillian into the common periodic-generator
representation. The Fourier decomposition is applied independently to the left and right
operator factors of every action and to its scalar coefficient, so periodic dependence in a
Hamiltonian, collapse operator, jump operator, rate, or any combination is supported.

Use [`liouvillian`](@ref) to construct the time-dependent map before calling
[`floquet_expansion`](@ref). The result is the same native
`PeriodicGenerator{Liouvillian}` as a manually assembled periodic Liouvillian.

See also [`PeriodicGenerator`](@ref), [`floquet_expansion`](@ref).
"""
function harmonics(L::Liouvillian, w::Symbolics.Num, t::Symbolics.Num)
  out = Dict{Int,Liouvillian}()
  for ((left, right), coefficient) in term_pairs(L)
    left_harmonics = harmonics(left, w, t)
    right_harmonics = harmonics(right, w, t)

    for (left_harmonic, left_component) in component_pairs(left_harmonics),
      (right_harmonic, right_component) in component_pairs(right_harmonics)

      foreach_harmonic_phase(coefficient, w, t) do coefficient_harmonic, phase_coefficient
        harmonic = left_harmonic + right_harmonic + coefficient_harmonic

        haskey(out, harmonic) || (out[harmonic] = zero(L))
        return add_term!(out[harmonic], left_component, right_component, phase_coefficient)
      end
    end
  end
  return PeriodicGenerator(out, w, zero(L))
end

function Base.show(io::IO, L::Liouvillian)
  return print(io, isempty(L) ? "0" : "Liouvillian($(length(L.terms)) terms)")
end

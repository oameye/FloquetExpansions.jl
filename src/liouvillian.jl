const LiouvillianAction = Tuple{SQA.QAdd,SQA.QAdd}
const LiouvillianTerms = Dict{LiouvillianAction,SQA.CNum}
const LiouvillianScalar = Union{Number,Symbolics.Num,SQA.CNum}
const LiouvillianChannelCollection = Union{Tuple,AbstractVector}

abstract type LiouvillianChannel end

struct CollapseChannel{O<:SQA.QField} <: LiouvillianChannel
  operator::O
end

struct RateWeightedJump{O<:SQA.QField,R<:LiouvillianScalar} <: LiouvillianChannel
  operator::O
  rate::R
end

"""
    Liouvillian

A symbolic linear map on density operators represented as a collected sum of elementary
actions ``ρ ↦ AρB``. The operator factors are SQA expressions and the scalar coefficients are
symbolic SQA coefficients.

Use [`hamiltonian_action`](@ref), [`dissipator`](@ref), or the channel keyword constructor to
build Liouvillians. Finite-order Floquet expansions can be more general than a generator in
explicit Lindblad form.
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

function action(left::SQA.QField, right::SQA.QField, coefficient::LiouvillianScalar=1)
  return raw_liouvillian(
    LiouvillianTerms((qadd(left), qadd(right)) => convert(SQA.CNum, coefficient))
  )
end

"""
    hamiltonian_action(H)

Construct the coherent density-operator action ``-i[H, ⋅]``.
"""
function hamiltonian_action(H::SQA.QField)
  Hq = qadd(H)
  identity = one(Hq)
  return action(Hq, identity, -im) + action(identity, Hq, im)
end

"""
    dissipator(L)

Construct the symbolic Lindblad dissipator
``D[L](ρ) = LρL† - (L†Lρ + ρL†L)/2``.
"""
function dissipator(L::SQA.QField)
  Lq = qadd(L)
  identity = one(Lq)
  norm = adjoint(Lq) * Lq
  return action(Lq, adjoint(Lq)) +
         action(norm, identity, -1 // 2) +
         action(identity, norm, -1 // 2)
end

"""
    Liouvillian(H; channels=())

Construct ``-i[H, ⋅] + Σ D[C] + Σ γD[J]``. Use [`collapse`](@ref) for complete collapse
operators and [`jump`](@ref) for paired bare jump operators and scalar rates. Rates are ordinary
symbolic coefficients; no positivity condition is inferred or certified.
"""
function Liouvillian(H::SQA.QField; channels::LiouvillianChannelCollection=())
  generator = hamiltonian_action(H)
  for channel in channels
    generator = generator + _channel_liouvillian(channel)
  end
  return generator
end

@inline _channel_liouvillian(channel::CollapseChannel) = dissipator(channel.operator)
@inline _channel_liouvillian(channel::RateWeightedJump) =
  channel.rate * dissipator(channel.operator)

"""
    collapse(operator::QField) -> CollapseChannel

Represent a complete collapse operator. It contributes ``D[operator]`` to a [`Liouvillian`](@ref)
with unit channel weight. Any amplitude or phase belonging to the collapse process is included in
`operator`.
"""
function collapse(operator::SQA.QField)
  return CollapseChannel(operator)
end

"""
    jump(operator::QField, rate) -> RateWeightedJump

Represent a bare jump operator with a separate symbolic rate. It contributes
``rate D[operator]`` to a [`Liouvillian`](@ref). The rate is not folded into the operator.
"""
function jump(operator::SQA.QField, rate::LiouvillianScalar)
  return RateWeightedJump(operator, rate)
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
    compose(A, B)

Compose Liouvillian maps so that `B` acts first and `A` acts second. For elementary actions
`(Aₗ, Aᵣ)` and `(Bₗ, Bᵣ)`, the resulting action is `(AₗBₗ, BᵣAᵣ)`.
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
    harmonics(L::Liouvillian, w, t) -> PeriodicGenerator{Liouvillian}

Lower a symbolic time-dependent Liouvillian into the common periodic-generator
representation. The Fourier decomposition is applied independently to the left and right
operator factors of every action and to its scalar coefficient, so periodic dependence in a
Hamiltonian, collapse operator, jump operator, rate, or any combination is supported.

Construct the time-dependent map with [`Liouvillian`](@ref), then pass it here before calling
[`floquet_expansion`](@ref). The result is the same native `PeriodicGenerator{Liouvillian}` as a
manually assembled periodic Liouvillian.
"""
function harmonics(L::Liouvillian, w::Symbolics.Num, t::Symbolics.Num)
  out = Dict{Int,Liouvillian}()
  for ((left, right), coefficient) in term_pairs(L)
    left_harmonics = harmonics(left, w, t)
    right_harmonics = harmonics(right, w, t)

    for (left_harmonic, left_component) in component_pairs(left_harmonics),
      (right_harmonic, right_component) in component_pairs(right_harmonics),
      phase_term in SQA.phase_terms(coefficient)

      phase_coefficient = phase_term.amplitude
      m, offset = harmonic_index(phase_term.phase, w, t)
      issymzero(offset) || (phase_coefficient *= SQA.expim(offset))
      harmonic = left_harmonic + right_harmonic + m

      haskey(out, harmonic) || (out[harmonic] = zero(L))
      add_term!(out[harmonic], left_component, right_component, phase_coefficient)
    end
  end
  return PeriodicGenerator(out, w, zero(L))
end

function Base.show(io::IO, L::Liouvillian)
  return print(io, isempty(L) ? "0" : "Liouvillian($(length(L.terms)) terms)")
end

const LiouvillianAction = Tuple{SQA.QAdd,SQA.QAdd}
const LiouvillianTerms = Dict{LiouvillianAction,SQA.CNum}
const LiouvillianScalar = Union{Number,Symbolics.Num,SQA.CNum}

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

function monomial_operator(term::SQA.QTerm)::SQA.QAdd
  arguments = SQA.QTermDict()
  arguments[term] = convert(SQA.CNum, 1)
  return SQA.QAdd(arguments, SQA.Index[])
end

@inline function add_action!(L::Liouvillian, key::LiouvillianAction, coefficient::SQA.CNum)
  updated = get(L.terms, key, convert(SQA.CNum, 0)) + coefficient
  iszero(updated) ? delete!(L.terms, key) : (L.terms[key] = updated)
  return L
end

function add_monomial_actions!(
  L::Liouvillian, left::SQA.QAdd, right::SQA.QAdd, coefficient::SQA.CNum
)
  for (left_term, left_coefficient) in SQA.expand_completeness(left),
    (right_term, right_coefficient) in SQA.expand_completeness(right)

    add_action!(
      L,
      (monomial_operator(left_term), monomial_operator(right_term)),
      coefficient * left_coefficient * right_coefficient,
    )
  end
  return L
end

@inline has_bound_sums(left::SQA.QAdd, right::SQA.QAdd) =
  !(isempty(left.indices) && isempty(right.indices))

function add_term!(L::Liouvillian, left::SQA.QAdd, right::SQA.QAdd, coefficient::SQA.CNum)
  (iszero(left) || iszero(right) || iszero(coefficient)) && return L
  has_bound_sums(left, right) && return add_action!(L, (left, right), coefficient)
  return add_monomial_actions!(L, left, right, coefficient)
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

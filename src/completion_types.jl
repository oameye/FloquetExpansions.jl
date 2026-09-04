"""
    Completion

Type-state hierarchy for the finite realization carried by a [`FloquetExpansion`](@ref).
Raw perturbative expansions use [`Uncompleted`](@ref); positive-completion algorithms replace
that state without rewriting the retained Floquet coefficients or micromotion.
"""
abstract type Completion end

"""
    Uncompleted <: Completion

Completion state of a raw Floquet expansion. Retained effective-generator and micromotion
components are the perturbative data produced by [`floquet_expansion`](@ref).
"""
struct Uncompleted <: Completion end

"""
    CompletionAlgorithm

Algorithm family for constructing a completely-positive finite realization of a
Liouvillian Floquet expansion.
"""
abstract type CompletionAlgorithm end

"""
    Gram <: CompletionAlgorithm

Algebraic Gram-factor completion algorithm.
"""
struct Gram <: CompletionAlgorithm end

"""
    Spectral <: CompletionAlgorithm

Spectral/HCM completion algorithm.
"""
struct Spectral <: CompletionAlgorithm end

"""
    CompletionFactorization

Common supertype for algorithm-specific positive-completion factorizations.
"""
abstract type CompletionFactorization end

# Internal microscopic provenance. These types deliberately live outside Liouvillian and
# PeriodicGenerator: generic map/Fourier algebra must not acquire physical channel metadata.
abstract type FloquetProvenance end
struct NoProvenance <: FloquetProvenance end

struct NonnegativeRateAssumption
  rate::SQA.CNum
end

@enum DissipativeSeedKind::UInt8 begin
  COLLAPSE_SEED = 0x01
  JUMP_SEED = 0x02
end

struct DissipativeSeedRef
  kind::DissipativeSeedKind
  index::Int
end

struct MicroscopicProvenance <: FloquetProvenance
  collapse_operators::Vector{SQA.QAdd}
  jump_operators::Vector{SQA.QAdd}
  jump_rates::Vector{SQA.CNum}
  rate_assumptions::Vector{NonnegativeRateAssumption}
  order::Vector{DissipativeSeedRef}
end

# Filled by the completion algorithms in later implementation stages. Every field is
# concretely parameterized so storing a completed state cannot introduce abstract fields.
struct PositiveCompletion{A<:CompletionAlgorithm,F,K,C,P,R,X<:CompletionFactorization,E} <:
       Completion
  algorithm::A
  frame::F
  kossakowski::K
  channels::C
  positivity_conditions::P
  regularity_conditions::R
  factorization::X
  generator::E
end

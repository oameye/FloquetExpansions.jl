"""
    Completion

Abstract completion state of a [`FloquetExpansion`](@ref). See [`positive_completion`](@ref).
"""
abstract type Completion end

"""
    Uncompleted <: Completion

State of a Floquet expansion before positive completion.
"""
struct Uncompleted <: Completion end

"""
    CompletionAlgorithm

Abstract algorithm selector for [`positive_completion`](@ref).
"""
abstract type CompletionAlgorithm end

"""
    Gram <: CompletionAlgorithm

Select algebraic Gram/Feshbach positive completion.

The construction determines a graded collapse-amplitude factor whose truncated Gram product
matches every retained Kossakowski coefficient. The untruncated finite Gram product then supplies
the positive continuation. Active dissipative directions are eliminated algebraically and the
construction recurses on the residual dark sector, so no Kossakowski eigendecomposition is
required.

`Gram()` can be used with any independent [`DissipativeFrame`](@ref). Supplying the frame
explicitly fixes the dissipative coordinates; otherwise [`positive_completion`](@ref) derives a
frame from the available Floquet data.
"""
struct Gram <: CompletionAlgorithm end

"""
    Spectral <: CompletionAlgorithm

Select perturbative spectral Haddadfarshi--Cui--Mintert (HCM) positive completion.

The construction follows perturbative decay-rate and branch-vector series in a dissipative frame
where the leading Kossakowski form is diagonal. Each admissible retained rate branch is completed
through a truncated square-root expansion, and the finite squared rates reconstruct a positive
Kossakowski form while preserving all retained coefficients. This is the completely-positive
high-frequency construction of Haddadfarshi, Cui, and Mintert [Haddadfarshi2015](@cite).

`Spectral()` therefore requires a suitable leading spectral frame. The two-argument
`positive_completion(expansion, Spectral())` uses the automatically derived frame only when that
frame already satisfies this condition; it does not diagonalize the leading Kossakowski form.
Degenerate leading sectors must first be resolved by an adapted degenerate perturbative basis
rather than by nondegenerate branch recursion.
"""
struct Spectral <: CompletionAlgorithm end

"""
    CompletionFactorization

Abstract type for factorization data produced by a positive-completion algorithm.
"""
abstract type CompletionFactorization end

"""
    CompletionObstruction

Raised when a retained dissipative leading form contains a direction incompatible with a
positive continuation on the current symbolic stratum.
"""
struct CompletionObstruction <: Exception
  rate_order::Int
  obstruction::SQA.CNum
  reason::Symbol
end

function Base.showerror(io::IO, error::CompletionObstruction)
  return print(
    io,
    "positive-completion obstruction at rate order ",
    error.rate_order,
    " (",
    error.reason,
    "): ",
    error.obstruction,
  )
end

"""
    FractionalJumpOnset

Raised when a positive dissipative rate first appears at odd inverse-frequency order, so an
integer-power collapse-amplitude representation would require a half-integer onset.
"""
struct FractionalJumpOnset <: Exception
  rate_order::Int
end

function Base.showerror(io::IO, error::FractionalJumpOnset)
  return print(
    io,
    "positive completion requires a fractional/Puiseux collapse-amplitude onset at rate order ",
    error.rate_order,
  )
end

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

struct PositiveCompletion{
  A<:CompletionAlgorithm,F,RK,K,C,P,R,X<:CompletionFactorization,H,E
} <: Completion
  algorithm::A
  frame::F
  retained_kossakowski::RK
  kossakowski::K
  channels::C
  positivity_conditions::P
  regularity_conditions::R
  factorization::X
  hamiltonian::H
  generator::E
end

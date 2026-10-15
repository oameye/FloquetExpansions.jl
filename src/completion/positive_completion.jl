function with_completion(expansion::FloquetExpansion, completion::PositiveCompletion)
  return FloquetExpansion(
    getfield(expansion, :generator),
    getfield(expansion, :micromotion_components),
    getfield(expansion, :effective_components),
    getfield(expansion, :gauge),
    getfield(expansion, :order),
    completion,
    getfield(expansion, :provenance),
  )
end

function finalize_positive_completion(
  expansion::FloquetExpansion,
  algorithm::CompletionAlgorithm,
  basis::JumpBasis,
  retained::RetainedGKSLData,
  completed_matrix::KossakowskiMatrix,
  completed_channels,
  conditions::CompletionConditions,
  factorization::CompletionFactorization,
)
  matrix_is_hermitian(completed_matrix) || throw(
    ArgumentError(
      "completed $(nameof(typeof(algorithm))) Kossakowski matrix is not Hermitian"
    ),
  )

  drive_frequency = getfield(expansion, :generator).wd
  retained_kossakowski = physical_kossakowski_series(retained.kossakowski, drive_frequency)
  coherent = retained.coherent
  lindbladian = SQA.liouvillian(
    coherent, first.(completed_channels); rates=last.(completed_channels)
  )
  generator = SQA.expand_completeness(SQA.SuperOperator(lindbladian))

  completion = PositiveCompletion(
    algorithm,
    basis,
    retained_kossakowski,
    completed_matrix,
    lindbladian,
    condition_coefficients(conditions.positivity),
    condition_coefficients(conditions.regularity),
    factorization,
    coherent,
    generator,
  )
  return with_completion(expansion, completion)
end

function positive_completion_impl(
  ::FloquetExpansion{G,P,E,Uncompleted,R}, algorithm::CompletionAlgorithm
) where {G<:Gauge,P<:PeriodicGenerator,E,R<:FloquetProvenance}
  return throw(
    ArgumentError("$(nameof(typeof(algorithm))) positive completion is not implemented yet")
  )
end

function positive_completion_impl(
  ::FloquetExpansion{G,P,E,Uncompleted,R}, algorithm::CompletionAlgorithm, ::JumpBasis
) where {G<:Gauge,P<:PeriodicGenerator,E,R<:FloquetProvenance}
  return throw(
    ArgumentError("$(nameof(typeof(algorithm))) positive completion is not implemented yet")
  )
end

function positive_completion_impl(
  expansion::FloquetExpansion{G,P,E,Uncompleted,R}, algorithm::Gram
) where {G<:Gauge,P<:PeriodicGenerator,E,R<:FloquetProvenance}
  return gram_positive_completion(expansion, automatic_jump_basis(expansion), algorithm)
end

function positive_completion_impl(
  expansion::FloquetExpansion{G,P,E,Uncompleted,R}, algorithm::Gram, basis::JumpBasis
) where {G<:Gauge,P<:PeriodicGenerator,E,R<:FloquetProvenance}
  return gram_positive_completion(expansion, basis, algorithm)
end

"""
    positive_completion(expansion::FloquetExpansion, algorithm::CompletionAlgorithm)
    positive_completion(expansion::FloquetExpansion, algorithm::CompletionAlgorithm,
                        basis::JumpBasis)

Return a positively completed Liouvillian Floquet expansion.

The returned [`FloquetExpansion`](@ref) preserves the retained effective components and
micromotion. [`effective_generator`](@ref) returns the completed finite generator, while
[`effective_component`](@ref) and [`micromotion`](@ref) continue to return the retained
high-frequency data.

The two-argument form derives a jump basis from the available Floquet data and, when
present, microscopic channel information. Pass a `JumpBasis` explicitly to fix the
operator representation, for example when comparing completion methods in the same physical basis.
Algorithm-specific basis requirements still apply; see [`Gram`](@ref) and [`Spectral`](@ref).
An expansion without dissipative terms completes in an empty derived basis with no channels; the
completed generator is then the Hamiltonian action of the retained coherent Hamiltonian.

Positive completion is defined only for Liouvillian expansions. Calling it on an already
completed expansion is an error. Use [`Gram`](@ref) or [`Spectral`](@ref) to select the algorithm.

# Examples

The following driven qubit first lets `Gram()` derive its jump basis, then uses an explicit
Cartesian basis for `Gram()` and a jump-adapted diagonal basis for `Spectral()`:

```jldoctest
julia> qubit = PauliSpace(:completion_doc);

julia> σx = Pauli(qubit, :σ, 1); σy = Pauli(qubit, :σ, 2); σz = Pauli(qubit, :σ, 3);

julia> σminus = (1 // 2) * (σx - im * σy);

julia> @variables ω::Real t::Real E::Real γ::Real;

julia> H = (1 // 2) * σz + E * cos(ω * t) * σx;

julia> raw = floquet_expansion(hamiltonian_action(H) + γ * dissipator(σminus), ω, t, VanVleck(), 3);

julia> cartesian = JumpBasis(σx, σy, σz);

julia> cp = positive_completion(raw, Gram(), cartesian)
FloquetExpansion{VanVleck{HoriDeprit}} of order 3

julia> kossakowski(cp)
3×3 Matrix{SecondQuantizedAlgebra.Coeff}:
 (1//4)*γ                                  …  0
 im*(((1//4)*(E^2)*γ) / (ω^2) - (1//4)*γ)     0
 0                                            ((1//2)*(E^2)*γ) / (ω^2)
```
"""
function positive_completion(
  expansion::FloquetExpansion{G,P,E,Uncompleted,R}, algorithm::CompletionAlgorithm
) where {G,P<:PeriodicGenerator{SQA.SuperOperator},E<:SQA.SuperOperator,R}
  return positive_completion_impl(expansion, algorithm)
end

function positive_completion(
  expansion::FloquetExpansion{G,P,E,Uncompleted,R},
  algorithm::CompletionAlgorithm,
  basis::JumpBasis,
) where {G,P<:PeriodicGenerator{SQA.SuperOperator},E<:SQA.SuperOperator,R}
  return positive_completion_impl(expansion, algorithm, basis)
end

function positive_completion(
  ::FloquetExpansion{G,P,E,Uncompleted,R}, ::CompletionAlgorithm
) where {G,P<:PeriodicGenerator{SQA.QAdd},E<:SQA.QAdd,R}
  return throw(
    ArgumentError("positive completion requires a Liouvillian Floquet expansion")
  )
end

function positive_completion(
  ::FloquetExpansion{G,P,E,Uncompleted,R}, ::CompletionAlgorithm, ::JumpBasis
) where {G,P<:PeriodicGenerator{SQA.QAdd},E<:SQA.QAdd,R}
  return throw(
    ArgumentError("positive completion requires a Liouvillian Floquet expansion")
  )
end

function positive_completion(
  ::FloquetExpansion{G,P,E,C,R}, ::CompletionAlgorithm
) where {G,P,E,C<:PositiveCompletion,R}
  return throw(ArgumentError("the Floquet expansion is already positively completed"))
end

function positive_completion(
  ::FloquetExpansion{G,P,E,C,R}, ::CompletionAlgorithm, ::JumpBasis
) where {G,P,E,C<:PositiveCompletion,R}
  return throw(ArgumentError("the Floquet expansion is already positively completed"))
end

function positive_completion(
  ::FloquetExpansion{G,P,E,C,R}, ::CompletionAlgorithm
) where {G,P,E,C<:NativeRealization,R}
  return throw(
    ArgumentError(
      "a GKSLNormalForm expansion is already of GKSL form; positive completion does not apply",
    ),
  )
end

function positive_completion(
  ::FloquetExpansion{G,P,E,C,R}, ::CompletionAlgorithm, ::JumpBasis
) where {G,P,E,C<:NativeRealization,R}
  return throw(
    ArgumentError(
      "a GKSLNormalForm expansion is already of GKSL form; positive completion does not apply",
    ),
  )
end

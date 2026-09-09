function with_completion(expansion::FloquetExpansion, completion::PositiveCompletion)
  return FloquetExpansion(
    getfield(expansion, :generator),
    getfield(expansion, :kick_components),
    getfield(expansion, :effective_components),
    getfield(expansion, :gauge),
    getfield(expansion, :order),
    completion,
    getfield(expansion, :provenance),
  )
end

@inline function stored_completion(expansion::FloquetExpansion)
  return getfield(expansion, :completion)
end

@inline stored_dissipative_frame(expansion::FloquetExpansion) =
  stored_completion(expansion).frame
@inline stored_retained_kossakowski(expansion::FloquetExpansion) =
  stored_completion(expansion).retained_kossakowski
@inline stored_completion_hamiltonian(expansion::FloquetExpansion) =
  stored_completion(expansion).hamiltonian

function finalize_positive_completion(
  expansion::FloquetExpansion,
  algorithm::CompletionAlgorithm,
  frame::DissipativeFrame,
  raw_matrices::Vector{KossakowskiMatrix},
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
  retained_kossakowski = physical_kossakowski_series(raw_matrices, drive_frequency)
  raw_generator = effective_generator(expansion)
  coherent = hamiltonian(raw_generator, frame)
  generator = liouvillian(coherent; channels=completed_channels)

  # The completed result owns its representation data. In particular, mutating a frame that
  # the caller supplied after completion must not invalidate cached coordinates or make the
  # stored generator inconsistent with the stored Kossakowski form.
  owned_frame = copy(frame)
  completion = PositiveCompletion(
    algorithm,
    owned_frame,
    retained_kossakowski,
    completed_matrix,
    completed_channels,
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
  ::FloquetExpansion{G,P,E,Uncompleted,R},
  algorithm::CompletionAlgorithm,
  ::DissipativeFrame,
) where {G<:Gauge,P<:PeriodicGenerator,E,R<:FloquetProvenance}
  return throw(
    ArgumentError("$(nameof(typeof(algorithm))) positive completion is not implemented yet")
  )
end

function positive_completion_impl(
  expansion::FloquetExpansion{G,P,E,Uncompleted,R}, algorithm::Gram
) where {G<:Gauge,P<:PeriodicGenerator,E,R<:FloquetProvenance}
  return gram_positive_completion(
    expansion, automatic_dissipative_frame(expansion), algorithm
  )
end

function positive_completion_impl(
  expansion::FloquetExpansion{G,P,E,Uncompleted,R}, algorithm::Gram, frame::DissipativeFrame
) where {G<:Gauge,P<:PeriodicGenerator,E,R<:FloquetProvenance}
  return gram_positive_completion(expansion, frame, algorithm)
end

"""
    positive_completion(expansion::FloquetExpansion, algorithm::CompletionAlgorithm)
    positive_completion(expansion::FloquetExpansion, algorithm::CompletionAlgorithm,
                        frame::DissipativeFrame)

Return a positively completed Liouvillian Floquet expansion.

The returned [`FloquetExpansion`](@ref) preserves the retained effective components and
micromotion. [`effective_generator`](@ref) returns the completed finite generator, while
[`effective_component`](@ref) and [`micromotion`](@ref) continue to return the retained
high-frequency data.

The two-argument form lets the completion algorithm determine a dissipative frame from the
available Floquet data and, when present, microscopic channel information. Pass a
[`DissipativeFrame`](@ref) explicitly to fix the representation. Automatic frame discovery is a
symbolic convenience frontend with runtime-dependent arity; the explicit-frame form is the
inference-oriented computational core.

Positive completion is defined only for Liouvillian expansions. Calling it on an already
completed expansion is an error. Use [`Gram`](@ref) or [`Spectral`](@ref) to select the algorithm.
"""
function positive_completion(
  expansion::FloquetExpansion{G,P,E,Uncompleted,R}, algorithm::CompletionAlgorithm
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian,R}
  return positive_completion_impl(expansion, algorithm)
end

function positive_completion(
  expansion::FloquetExpansion{G,P,E,Uncompleted,R},
  algorithm::CompletionAlgorithm,
  frame::DissipativeFrame,
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian,R}
  return positive_completion_impl(expansion, algorithm, frame)
end

function positive_completion(
  ::FloquetExpansion{G,P,E,Uncompleted,R}, ::CompletionAlgorithm
) where {G,P<:PeriodicGenerator{SQA.QAdd},E<:SQA.QAdd,R}
  return throw(
    ArgumentError("positive completion requires a Liouvillian Floquet expansion")
  )
end

function positive_completion(
  ::FloquetExpansion{G,P,E,Uncompleted,R}, ::CompletionAlgorithm, ::DissipativeFrame
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
  ::FloquetExpansion{G,P,E,C,R}, ::CompletionAlgorithm, ::DissipativeFrame
) where {G,P,E,C<:PositiveCompletion,R}
  return throw(ArgumentError("the Floquet expansion is already positively completed"))
end

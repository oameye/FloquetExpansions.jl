function _with_completion(expansion::FloquetExpansion, completion::PositiveCompletion)
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

function _positive_completion(::FloquetExpansion, algorithm::CompletionAlgorithm)
  return throw(
    ArgumentError("$(nameof(typeof(algorithm))) positive completion is not implemented yet")
  )
end

function _positive_completion(
  ::FloquetExpansion, algorithm::CompletionAlgorithm, ::DissipativeFrame
)
  return throw(
    ArgumentError("$(nameof(typeof(algorithm))) positive completion is not implemented yet")
  )
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

The two-argument form is available for expansions constructed from
`floquet_expansion(H, ωd, t, gauge, order; channels=...)`. For an expansion constructed from an
explicit [`Liouvillian`](@ref), [`harmonics`](@ref), or [`PeriodicGenerator`](@ref), pass a
[`DissipativeFrame`](@ref) explicitly.

Positive completion is defined only for Liouvillian expansions. Calling it on an already
completed expansion is an error. Use [`Gram`](@ref) or [`Spectral`](@ref) to select the algorithm.
"""
function positive_completion(
  expansion::FloquetExpansion{G,P,E,Uncompleted,MicroscopicProvenance},
  algorithm::CompletionAlgorithm,
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian}
  return _positive_completion(expansion, algorithm)
end

function positive_completion(
  ::FloquetExpansion{G,P,E,Uncompleted,NoProvenance}, ::CompletionAlgorithm
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian}
  return throw(
    ArgumentError(
      "automatic positive completion requires an expansion constructed with `channels`; pass a `DissipativeFrame` for an explicitly constructed Liouvillian",
    ),
  )
end

function positive_completion(
  expansion::FloquetExpansion{G,P,E,Uncompleted,R},
  algorithm::CompletionAlgorithm,
  frame::DissipativeFrame,
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian,R}
  return _positive_completion(expansion, algorithm, frame)
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

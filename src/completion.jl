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

Construct a completely-positive finite realization of a Liouvillian Floquet expansion.
The returned object is a new [`FloquetExpansion`](@ref): retained effective components,
micromotion, gauge, expansion order, and the input periodic generator are preserved.
Only the finite no-index [`effective_generator`](@ref) realization changes.

The overload without a frame uses microscopic channel provenance retained by the high-level
Hamiltonian-plus-channels constructor. A manually lowered Liouvillian has no such provenance and
therefore requires the explicit `frame` overload. Calling `positive_completion` on a Hamiltonian
expansion or on an already completed expansion is an error.

`Gram()` and `Spectral()` establish the public algorithm selection API; their mathematical
implementations are provided by the corresponding completion stages.
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
      "automatic positive completion requires microscopic channel provenance; pass a `DissipativeFrame` explicitly for a manually lowered Liouvillian",
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

@inline function stored_completion(expansion::FloquetExpansion)
  return getfield(expansion, :completion)
end

@inline stored_dissipative_frame(expansion::FloquetExpansion) =
  stored_completion(expansion).frame
@inline stored_retained_kossakowski(expansion::FloquetExpansion) =
  stored_completion(expansion).retained_kossakowski
@inline stored_completion_hamiltonian(expansion::FloquetExpansion) =
  stored_completion(expansion).hamiltonian

function effective_generator(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:Union{PositiveCompletion,NativeRealization},R}
  return getfield(expansion, :completion).generator
end

function hamiltonian(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {
  G,
  P<:PeriodicGenerator{Liouvillian},
  E<:Liouvillian,
  C<:Union{PositiveCompletion,NativeRealization},
  R<:FloquetProvenance,
}
  return stored_completion_hamiltonian(expansion)::SQA.QAdd
end

"""
    dissipative_frame(expansion::FloquetExpansion)

Return an independent copy of the ordered dissipative frame owned by a positively completed
Floquet expansion.
"""
function dissipative_frame(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:Union{PositiveCompletion,NativeRealization},R}
  return copy(stored_dissipative_frame(expansion))
end

"""
    channels(expansion::FloquetExpansion)

Return the completed physical channels. Together with [`hamiltonian`](@ref), these channels
reconstruct [`effective_generator`](@ref) exactly.
"""
function channels(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:Union{PositiveCompletion,NativeRealization},R}
  return copy(stored_completion(expansion).channels)
end

"""
    positivity_conditions(expansion::FloquetExpansion)

Return scalar conditions assumed nonnegative on the symbolic stratum used by positive
completion.
"""
function positivity_conditions(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:Union{PositiveCompletion,NativeRealization},R}
  return copy(stored_completion(expansion).positivity_conditions)
end

"""
    regularity_conditions(expansion::FloquetExpansion)

Return scalar conditions assumed nonzero to remain on the fixed-rank symbolic stratum used by
positive completion.
"""
function regularity_conditions(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:Union{PositiveCompletion,NativeRealization},R}
  return copy(stored_completion(expansion).regularity_conditions)
end

"""
    factorization(expansion::FloquetExpansion)

Return the algorithm-specific factorization data stored by a positively completed expansion.
For [`Gram`](@ref) completion this is a [`GramFactorization`](@ref).
"""
function factorization(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:PositiveCompletion,R}
  return copy(stored_completion(expansion).factorization)
end

function factorization(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:NativeRealization,R}
  return stored_completion(expansion).factorization
end

"""
    kossakowski(expansion::FloquetExpansion)

Return the finite completed Kossakowski matrix in [`dissipative_frame`](@ref). This no-frame
form is defined for positively completed Liouvillian expansions.
"""
function kossakowski(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {
  G,
  P<:PeriodicGenerator{Liouvillian},
  E<:Liouvillian,
  C<:Union{PositiveCompletion,NativeRealization},
  R,
}
  return copy(stored_completion(expansion).kossakowski)
end

"""
    kossakowski_component(expansion::FloquetExpansion, n::Int)

Return the cached retained order-`n` Kossakowski contribution in the completed expansion's
stored dissipative frame. Positive completion does not alter retained components.
"""
function kossakowski_component(
  expansion::FloquetExpansion{G,P,E,C,R}, n::Int
) where {
  G,
  P<:PeriodicGenerator{Liouvillian},
  E<:Liouvillian,
  C<:Union{PositiveCompletion,NativeRealization},
  R,
}
  retained = stored_retained_kossakowski(expansion)
  0 <= n < length(retained) || throw(BoundsError(retained, n + 1))
  return copy(retained[n + 1])
end

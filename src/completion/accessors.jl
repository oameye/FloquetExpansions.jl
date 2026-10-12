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

function SQA.hamiltonian(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {
  G,
  P<:PeriodicGenerator{SQA.SuperOperator},
  E<:SQA.SuperOperator,
  C<:Union{PositiveCompletion,NativeRealization},
  R<:FloquetProvenance,
}
  return stored_completion_hamiltonian(expansion)
end

"""
    dissipative_frame(expansion::FloquetExpansion)

Return the ordered dissipative frame of a positively completed or [`GKSLNormalForm`](@ref)
Floquet expansion.
"""
function dissipative_frame(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:Union{PositiveCompletion,NativeRealization},R}
  return stored_dissipative_frame(expansion)
end

"""
    Lindbladian(expansion::FloquetExpansion)

Return the effective generator of a positively completed or [`GKSLNormalForm`](@ref) Floquet
expansion as a SecondQuantizedAlgebra `Lindbladian`. Its rates are nonnegative on the stratum
of [`positivity_conditions`](@ref), and its `SuperOperator` lowering is
[`effective_generator`](@ref).

See also [`channels`](@ref), [`hamiltonian`](@ref).
"""
function SQA.Lindbladian(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:Union{PositiveCompletion,NativeRealization},R}
  return stored_completion(expansion).lindbladian
end

"""
    channels(expansion::FloquetExpansion)

Return the channels of [`Lindbladian`](@ref)`(expansion)` for a positively completed or
[`GKSLNormalForm`](@ref) Floquet expansion, as a vector of `channel => rate` pairs in the
deterministic order of `SecondQuantizedAlgebra.channel_order_key`. Read the jump operator of a
channel with `SecondQuantizedAlgebra.jump`.
"""
function SQA.channels(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:Union{PositiveCompletion,NativeRealization},R}
  return SQA.channels(stored_completion(expansion).lindbladian)
end

"""
    positivity_conditions(expansion::FloquetExpansion)

Return scalar conditions assumed nonnegative on the symbolic stratum used by positive
completion. A [`GKSLNormalForm`](@ref) expansion with symbolic parameters returns the signed
parameters and the drive frequency, and with numeric coefficients none.
"""
function positivity_conditions(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:Union{PositiveCompletion,NativeRealization},R}
  return copy(stored_completion(expansion).positivity_conditions)
end

"""
    regularity_conditions(expansion::FloquetExpansion)

Return scalar conditions assumed nonzero to remain on the fixed-rank symbolic stratum used by
positive completion. A [`GKSLNormalForm`](@ref) expansion with symbolic parameters returns the
irreducible denominators of its coefficients, and with numeric coefficients none.
"""
function regularity_conditions(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:Union{PositiveCompletion,NativeRealization},R}
  return copy(stored_completion(expansion).regularity_conditions)
end

"""
    factorization(expansion::FloquetExpansion)

Return the algorithm-specific factorization data stored by a positively completed expansion.
For [`Gram`](@ref) completion this is a [`GramFactorization`](@ref). A [`GKSLNormalForm`](@ref)
expansion returns internal expansion data, which is not part of the public interface.
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
form is defined for positively completed and [`GKSLNormalForm`](@ref) Liouvillian expansions.
"""
function SQA.kossakowski(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {
  G,
  P<:PeriodicGenerator{SQA.SuperOperator},
  E<:SQA.SuperOperator,
  C<:Union{PositiveCompletion,NativeRealization},
  R,
}
  return copy(stored_completion(expansion).kossakowski)
end

"""
    kossakowski_component(expansion::FloquetExpansion, n::Int)

Return the cached retained order-`n` Kossakowski contribution in the stored dissipative frame of
a positively completed or [`GKSLNormalForm`](@ref) expansion. Positive completion does not alter
retained components.
"""
function kossakowski_component(
  expansion::FloquetExpansion{G,P,E,C,R}, n::Int
) where {
  G,
  P<:PeriodicGenerator{SQA.SuperOperator},
  E<:SQA.SuperOperator,
  C<:Union{PositiveCompletion,NativeRealization},
  R,
}
  retained = stored_retained_kossakowski(expansion)
  0 <= n < length(retained) || throw(BoundsError(retained, n + 1))
  return copy(retained[n + 1])
end

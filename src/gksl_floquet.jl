"""
    kossakowski(expansion::FloquetExpansion, frame::DissipativeFrame)

Return the Hermitian Kossakowski matrix of the finite effective Liouvillian in the ordered
`frame`.

See also [`DissipativeFrame`](@ref), [`kossakowski_component`](@ref), [`hamiltonian`](@ref).
"""
function kossakowski(
  expansion::FloquetExpansion{G,P,E}, frame::DissipativeFrame
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian}
  generator = effective_generator(expansion)::E
  return kossakowski(generator, frame)::KossakowskiMatrix
end

"""
    kossakowski_component(expansion::FloquetExpansion, frame::DissipativeFrame, n::Int)

Return the order-`n` Kossakowski contribution of a Liouvillian Floquet expansion in the
ordered dissipative `frame`, including the corresponding inverse-drive-frequency scaling.

See also [`kossakowski`](@ref), [`effective_component`](@ref).
"""
function kossakowski_component(
  expansion::FloquetExpansion{G,P,E}, frame::DissipativeFrame, n::Int
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian}
  generator = effective_component(expansion, n)::E
  return kossakowski(generator, frame)::KossakowskiMatrix
end

"""
    hamiltonian(expansion::FloquetExpansion)

Return the finite effective Hamiltonian of a Hamiltonian Floquet expansion, or the coherent
Hamiltonian sector of a Liouvillian Floquet expansion. The latter is defined modulo an
additive multiple of the identity.

See also [`hamiltonian_component`](@ref), [`effective_generator`](@ref).
"""
function hamiltonian(
  expansion::FloquetExpansion{G,P,E}
) where {G,P<:PeriodicGenerator{SQA.QAdd},E<:SQA.QAdd}
  return effective_generator(expansion)::E
end

function hamiltonian(
  expansion::FloquetExpansion{G,P,E}
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian}
  generator = effective_generator(expansion)::E
  return hamiltonian(generator)::SQA.QAdd
end

"""
    hamiltonian_component(expansion::FloquetExpansion, n::Int)

Return the coherent order-`n` contribution of a Floquet expansion, including its
inverse-drive-frequency scaling. The Hamiltonian is defined modulo an additive scalar
multiple of the identity.

See also [`hamiltonian`](@ref), [`effective_component`](@ref).
"""
function hamiltonian_component(
  expansion::FloquetExpansion{G,P,E}, n::Int
) where {G,P<:PeriodicGenerator{SQA.QAdd},E<:SQA.QAdd}
  return effective_component(expansion, n)::E
end

function hamiltonian_component(
  expansion::FloquetExpansion{G,P,E}, n::Int
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian}
  generator = effective_component(expansion, n)::E
  return hamiltonian(generator)::SQA.QAdd
end

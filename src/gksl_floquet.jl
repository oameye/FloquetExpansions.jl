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

# Examples

Construct a frame before extracting the finite form and one retained component:

```jldoctest
julia> qubit = PauliSpace(:kossakowski_doc);

julia> σx = Pauli(qubit, :σ, 1); σy = Pauli(qubit, :σ, 2); σz = Pauli(qubit, :σ, 3);

julia> σminus = (1 // 2) * (σx - im * σy);

julia> @variables ω::Real t::Real E::Real γ::Real;

julia> H = (1 // 2) * σz + E * cos(ω * t) * σx;

julia> raw = floquet_expansion(H, ω, t, VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))), 3; channels=(jump(σminus, γ),));

julia> frame = DissipativeFrame(σx, σy, σz);

julia> d = kossakowski_component(raw, frame, 1)
3×3 Matrix{SecondQuantizedAlgebra.Coeff}:
 0  0  0
 0  0  0
 0  0  0

julia> d = kossakowski_component(raw, frame, 2)
3×3 Matrix{SecondQuantizedAlgebra.Coeff}:
 0                            …  0
 ((1//4)*(E^2)*im*γ) / (ω^2)     0
 0                               ((1//2)*(E^2)*γ) / (ω^2)
```

See also [`DissipativeFrame`](@ref), [`kossakowski`](@ref), [`effective_component`](@ref).
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

For a positively completed Liouvillian expansion, this coherent sector is unchanged from the
retained Floquet expansion; positive completion modifies only dissipative information beyond the
retained order.

See also [`hamiltonian_component`](@ref), [`effective_generator`](@ref).
"""
function hamiltonian(
  expansion::FloquetExpansion{G,P,E}
) where {G,P<:PeriodicGenerator{SQA.QAdd},E<:SQA.QAdd}
  return effective_generator(expansion)::E
end

function hamiltonian(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {
  G,
  P<:PeriodicGenerator{Liouvillian},
  E<:Liouvillian,
  C<:PositiveCompletion,
  R<:FloquetProvenance,
}
  return stored_completion_hamiltonian(expansion)::SQA.QAdd
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

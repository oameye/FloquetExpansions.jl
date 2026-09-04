function kossakowski(
  expansion::FloquetExpansion{G,P,E}, frame::DissipativeFrame
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian}
  generator = effective_generator(expansion)::E
  return kossakowski(generator, frame)::KossakowskiMatrix
end

function kossakowski_component(
  expansion::FloquetExpansion{G,P,E}, frame::DissipativeFrame, n::Int
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian}
  generator = effective_generator(expansion, n)::E
  return kossakowski(generator, frame)::KossakowskiMatrix
end

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

function hamiltonian_component(
  expansion::FloquetExpansion{G,P,E}, n::Int
) where {G,P<:PeriodicGenerator{SQA.QAdd},E<:SQA.QAdd}
  return effective_generator(expansion, n)::E
end

function hamiltonian_component(
  expansion::FloquetExpansion{G,P,E}, n::Int
) where {G,P<:PeriodicGenerator{Liouvillian},E<:Liouvillian}
  generator = effective_generator(expansion, n)::E
  return hamiltonian(generator)::SQA.QAdd
end

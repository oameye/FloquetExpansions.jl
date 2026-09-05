using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

@variables ω::Real t::Real Ω::Real

@testset "transverse drive opens a second-order dephasing channel" begin
  pauli = PauliSpace(:gram_recursive_pauli)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σx, σy, σz)

  H = Ω * cos(ω * t) * σx
  expansion = floquet_expansion(H, ω, t, VanVleck(), 3; channels=(collapse(σz),))
  completion = positive_completion(expansion, Gram(), frame)
  gram = factorization(completion)

  @test gram.onsets == [0, 1]
  @test length(gram.stages) >= 2
  @test gram.stages[1].grade == 0
  @test gram.stages[2].grade == 2
  @test size(gram.amplitudes[1], 2) == 2
  @test all(iszero, gram.amplitudes[1][:, 2])
  @test any(!iszero(value) for value in gram.amplitudes[2][:, 2])
  @test liouvillian(hamiltonian(completion); channels=channels(completion)) ==
    effective_generator(completion)
end

using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

function completion_matrix_equal(left, right)
  size(left) == size(right) || return false
  return all(iszero(SQA.simplify(left[index] - right[index])) for index in eachindex(left))
end

@testset "shared CP-completion physical invariants" begin
  pauli = PauliSpace(:cp_validation_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  @variables ω::Real t::Real Ω::Real

  expansion = floquet_expansion(
    Ω * cos(ω * t) * σx, ω, t, VanVleck(), 3; channels=(collapse(σz),)
  )

  for method in (Gram(), Spectral())
    completion = positive_completion(expansion, method, frame)

    for n in 0:2
      @test effective_component(completion, n) == effective_component(expansion, n)
      @test hamiltonian_component(completion, n) == hamiltonian_component(expansion, n)
      @test completion_matrix_equal(
        kossakowski_component(completion, n), kossakowski_component(expansion, frame, n)
      )
    end

    @test micromotion(completion) == micromotion(expansion)
    @test SQA.simplify(liouvillian(hamiltonian(completion); channels=channels(completion))) ==
      SQA.simplify(effective_generator(completion))
  end
end

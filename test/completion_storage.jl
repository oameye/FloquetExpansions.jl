using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

@testset "completed representation owns frame and cached retained data" begin
  fock = FockSpace(:completion_storage)
  a = Destroy(fock, :a)
  @variables ω::Real t::Real

  frame = DissipativeFrame(a)
  expansion = floquet_expansion(0 * a, ω, t, VanVleck(), 1; channels=(collapse(a),))
  completion = @inferred positive_completion(expansion, Gram(), frame)

  retained = kossakowski_component(completion, 0)
  completed = kossakowski(completion)
  stored_frame = dissipative_frame(completion)

  # Mutating the caller-owned frame after completion must not alter the finalized result.
  frame.coordinates[1, 1] = convert(SQA.CNum, 0)
  frame.pivot_inverse[1, 1] = convert(SQA.CNum, 0)
  empty!(frame.pivot_rows)

  @test dissipative_frame(completion) == stored_frame
  @test kossakowski_component(completion, 0) == retained
  @test kossakowski(completion) == completed
  @test liouvillian(hamiltonian(completion); channels=channels(completion)) ==
    effective_generator(completion)

  # Public access remains defensive while internal accessors use the owned representation.
  copied = dissipative_frame(completion)
  copied.coordinates[1, 1] = convert(SQA.CNum, 0)
  @test !iszero(dissipative_frame(completion).coordinates[1, 1])

  cached = kossakowski_component(completion, 0)
  cached[1, 1] = convert(SQA.CNum, 0)
  @test !iszero(kossakowski_component(completion, 0)[1, 1])
end

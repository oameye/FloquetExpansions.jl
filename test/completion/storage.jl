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
  expansion = floquet_expansion(dissipator(a), ω, t, VanVleck(), 1)
  completion = @inferred positive_completion(expansion, Gram(), frame)

  retained = kossakowski_component(completion, 0)
  completed = kossakowski(completion)

  @test dissipative_frame(completion) == frame
  @test SQA.expand_completeness(SQA.SuperOperator(Lindbladian(completion))) ==
    effective_generator(completion)

  retained[1, 1] = convert(SQA.CNum, 0)
  completed[1, 1] = convert(SQA.CNum, 0)
  @test !iszero(kossakowski_component(completion, 0)[1, 1])
  @test !iszero(kossakowski(completion)[1, 1])
end

using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

h = FockSpace(:completion_state)
a = Destroy(h, :a)
@variables ω::Real t::Real γ::Real

@testset "completion diagnostics have semantic displays" begin
  expansion = floquet_expansion(γ * dissipator(a), ω, t, VanVleck(), 1)
  frame = JumpBasis(a)
  gram = factorization(positive_completion(expansion, Gram(), frame))
  spectral = factorization(positive_completion(expansion, Spectral(), frame))

  @test sprint(show, gram) == "GramFactorization(1 collapse channel, onsets=[0], stages=1)"
  @test sprint(show, MIME"text/plain"(), gram) ==
    "GramFactorization with 1 collapse channel\n" *
        "  onsets: [0]\n" *
        "  stages: 1\n" *
        "    GramStage(grade=0, active_rank=1, dark_rank=0)"
  @test sprint(show, spectral) ==
    "SpectralFactorization(1 branch, onsets=[0], puiseux=[false])"
  @test sprint(show, MIME"text/plain"(), spectral) ==
    "SpectralFactorization with 1 branch\n  onsets: [0]\n  puiseux: [false]"
end

@testset "raw Floquet expansions carry uncompleted state" begin
  H = a' * a + cos(ω * t) * (a + a')
  vv = @inferred floquet_expansion(H, ω, t, VanVleck(), 2)

  @test vv.completion isa Uncompleted
  @test @inferred(effective_component(vv, 0)) == hamiltonian_component(vv, 0)
  @test_throws ArgumentError effective_component(vv, -1)
  @test_throws ArgumentError effective_component(vv, vv.order)
  @test_throws MethodError effective_generator(vv, 0)
end

@testset "physical and lowered open-system constructors agree" begin
  H = a' * a + cos(ω * t) * (a + a')
  L = hamiltonian_action(H) + dissipator(a) + γ * dissipator(a')

  physical = @inferred floquet_expansion(L, ω, t, VanVleck(), 2)
  lowered = @inferred floquet_expansion(SQA.SuperOperator(L), ω, t, VanVleck(), 2)

  @test physical.completion isa Uncompleted
  @test lowered.completion isa Uncompleted
  @test effective_generator(physical) == effective_generator(lowered)
  @test effective_component(physical, 0) == effective_component(lowered, 0)
  @test effective_component(physical, 1) == effective_component(lowered, 1)
  @test micromotion(physical) == micromotion(lowered)
end

@testset "Lindbladian input requires a Hermitian Hamiltonian" begin
  @test_throws ArgumentError floquet_expansion(
    hamiltonian_action(a) + dissipator(a), ω, t, VanVleck(), 1
  )
end

@testset "positive-completion algorithms establish the public dispatch boundary" begin
  H = a' * a + cos(ω * t) * (a + a')
  physical_L = hamiltonian_action(H) + γ * dissipator(a)
  physical = floquet_expansion(physical_L, ω, t, VanVleck(), 1)
  lowered = @inferred floquet_expansion(SQA.SuperOperator(physical_L), ω, t, VanVleck(), 1)
  periodic = @inferred floquet_expansion(harmonics(physical_L, ω, t), VanVleck(), 1)
  coherent = floquet_expansion(H, ω, t, VanVleck(), 1)
  frame = JumpBasis(a)

  @test Gram() isa CompletionAlgorithm
  @test Spectral() isa CompletionAlgorithm
  @test applicable(positive_completion, physical, Gram())
  @test applicable(positive_completion, lowered, Gram())
  @test applicable(positive_completion, periodic, Gram())
  @test applicable(positive_completion, lowered, Spectral(), frame)
  @test applicable(positive_completion, periodic, Gram(), frame)

  @test_throws ArgumentError positive_completion(coherent, Gram())
  @test_throws ArgumentError positive_completion(coherent, Gram(), frame)
end

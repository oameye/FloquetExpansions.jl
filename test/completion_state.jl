using Test
using FloquetExpansions
using Symbolics: @variables

h = FockSpace(:completion_state)
a = Destroy(h, :a)
@variables ω::Real t::Real γ::Real

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
  channels = (collapse(a), jump(a', γ))

  physical = @inferred floquet_expansion(H, ω, t, VanVleck(), 2; channels)
  lowered = @inferred floquet_expansion(liouvillian(H; channels), ω, t, VanVleck(), 2)

  @test physical.completion isa Uncompleted
  @test lowered.completion isa Uncompleted
  @test effective_generator(physical) == effective_generator(lowered)
  @test effective_component(physical, 0) == effective_component(lowered, 0)
  @test effective_component(physical, 1) == effective_component(lowered, 1)
  @test micromotion(physical) == micromotion(lowered)
end

@testset "channel collection representation does not change raw dynamics" begin
  H = a' * a + cos(ω * t) * (a + a')
  tuple_channels = (collapse(a), jump(a', γ))
  vector_channels = [collapse(a), jump(a', γ)]

  tuple_vv = @inferred floquet_expansion(H, ω, t, VanVleck(), 2; channels=tuple_channels)
  vector_vv = @inferred floquet_expansion(H, ω, t, VanVleck(), 2; channels=vector_channels)
  reversed_vv = @inferred floquet_expansion(
    H, ω, t, VanVleck(), 2; channels=(jump(a', γ), collapse(a))
  )

  for other in (vector_vv, reversed_vv)
    @test effective_generator(other) == effective_generator(tuple_vv)
    @test effective_component(other, 0) == effective_component(tuple_vv, 0)
    @test effective_component(other, 1) == effective_component(tuple_vv, 1)
    @test micromotion(other) == micromotion(tuple_vv)
  end
end

@testset "physical Hamiltonian input remains a Hamiltonian" begin
  H = a' * a + cos(ω * t) * (a + a')
  empty_channels = [collapse(a)]
  empty!(empty_channels)

  @test_throws ArgumentError floquet_expansion(
    H, ω, t, VanVleck(), 1; channels=empty_channels
  )
  @test_throws ArgumentError floquet_expansion(
    a, ω, t, VanVleck(), 1; channels=(collapse(a),)
  )
end

@testset "positive-completion algorithms establish the public dispatch boundary" begin
  H = a' * a + cos(ω * t) * (a + a')
  physical = floquet_expansion(H, ω, t, VanVleck(), 1; channels=(jump(a, γ),))
  lowered_L = liouvillian(H; channels=(jump(a, γ),))
  lowered = @inferred floquet_expansion(lowered_L, ω, t, VanVleck(), 1)
  periodic = @inferred floquet_expansion(harmonics(lowered_L, ω, t), VanVleck(), 1)
  coherent = floquet_expansion(H, ω, t, VanVleck(), 1)
  frame = DissipativeFrame(a)

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

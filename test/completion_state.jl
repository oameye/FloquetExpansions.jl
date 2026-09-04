using Test
using FloquetExpansions
using Symbolics: @variables

h = FockSpace(:completion_state)
a = Destroy(h, :a)
@variables ω::Real t::Real γ::Real

function error_text(f)
  try
    f()
  catch err
    return sprint(showerror, err)
  end
  return ""
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

@testset "high-level dissipative construction remains concrete" begin
  H = a' * a + cos(ω * t) * (a + a')
  channels = [collapse(a), jump(a', γ)]
  vv = @inferred floquet_expansion(H, ω, t, VanVleck(), 2; channels)

  @test vv.completion isa Uncompleted
  @test @inferred(effective_generator(vv)) isa Liouvillian
  @test @inferred(effective_component(vv, 0)) isa Liouvillian
  @test @inferred(micromotion(vv)) isa PeriodicGenerator{Liouvillian}
end

@testset "positive-completion algorithms establish the public dispatch boundary" begin
  H = a' * a + cos(ω * t) * (a + a')
  physical = floquet_expansion(H, ω, t, VanVleck(), 1; channels=(jump(a, γ),))
  lowered_L = liouvillian(H; channels=(jump(a, γ),))
  lowered = @inferred floquet_expansion(lowered_L, ω, t, VanVleck(), 1)
  periodic = @inferred floquet_expansion(harmonics(lowered_L, ω, t), VanVleck(), 1)
  coherent = floquet_expansion(H, ω, t, VanVleck(), 1)
  frame = DissipativeFrame(a)

  @test lowered.completion isa Uncompleted
  @test periodic.completion isa Uncompleted
  @test Gram() isa CompletionAlgorithm
  @test Spectral() isa CompletionAlgorithm

  # The high-level physical path reaches the automatic-completion algorithm seam.
  @test occursin("not implemented", error_text(() -> positive_completion(physical, Gram())))

  # Explicitly constructed Liouvillians need a dissipative frame.
  @test occursin("DissipativeFrame", error_text(() -> positive_completion(lowered, Gram())))
  @test occursin(
    "not implemented", error_text(() -> positive_completion(lowered, Spectral(), frame))
  )

  # Positive completion is an open-system operation.
  @test occursin("Liouvillian", error_text(() -> positive_completion(coherent, Gram())))
  @test occursin(
    "Liouvillian", error_text(() -> positive_completion(coherent, Gram(), frame))
  )
end

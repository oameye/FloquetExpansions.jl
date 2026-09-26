using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const FE = FloquetExpansions
const SQA = SecondQuantizedAlgebra

space = NLevelSpace(:native_gksl_inference, 2)
σx = Transition(space, :σ, 1, 2) + Transition(space, :σ, 2, 1)
σz = Transition(space, :σ, 1, 1) - Transition(space, :σ, 2, 2)
σminus = Transition(space, :σ, 1, 2)
@variables w_native_inference::Real t_native_inference::Real γ_native_inference::Real

w = w_native_inference
t = t_native_inference
γ = γ_native_inference
H = σz + σx * SQA.expim(-w * t) + σx * SQA.expim(w * t)

@testset "native GKSL expansion amplitude core is inferred" begin
  seed = only(FE.jump_amplitude_seeds((jump(σminus, γ),), w, t))
  coherent = floquet_expansion(harmonics(H, w, t), VanVleck(), 2)
  K = coherent.micromotion_components

  transported = @inferred FE.transport_jump_amplitude(seed, K, 2)
  @test transported isa FE.TransportedJumpAmplitude

  channels = @inferred FE.harmonic_jump_channels([transported])
  @test channels isa Vector{FE.HarmonicJumpChannel}

  generator = @inferred FE.native_gksl_generator(coherent, channels)
  @test generator isa Liouvillian
end

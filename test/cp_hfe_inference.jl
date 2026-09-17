using Test
using FloquetExpansions
using JET: JET
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const FEInference = FloquetExpansions
const SQAInference = SecondQuantizedAlgebra

space_inference = NLevelSpace(:cp_hfe_inference, 2)
σx_inference = Transition(space_inference, :σ, 1, 2) + Transition(space_inference, :σ, 2, 1)
σz_inference = Transition(space_inference, :σ, 1, 1) - Transition(space_inference, :σ, 2, 2)
σminus_inference = Transition(space_inference, :σ, 1, 2)
@variables w_cp_inference::Real t_cp_inference::Real γ_cp_inference::Real

H_cp_inference =
  σz_inference +
  σx_inference * SQAInference.expim(-w_cp_inference * t_cp_inference) +
  σx_inference * SQAInference.expim(w_cp_inference * t_cp_inference)

@testset "native CP-HFE explicit amplitude core is inference-friendly" begin
  seed = only(
    FEInference.physical_amplitude_seeds(
      (jump(σminus_inference, γ_cp_inference),), w_cp_inference, t_cp_inference
    ),
  )
  coherent = floquet_expansion(
    harmonics(H_cp_inference, w_cp_inference, t_cp_inference), VanVleck(), 2
  )
  kicks = getfield(coherent, :kick_components)

  transported = @inferred FEInference.transport_amplitude_series(seed, kicks, 2)
  @test transported isa FEInference.TransportedAmplitudeSeries
  JET.@test_opt target_modules=(FloquetExpansions,) FEInference.transport_amplitude_series(
    seed, kicks, 2
  )

  rows = @inferred FEInference.reconstruct_cp_amplitude_channels([transported])
  @test rows isa Vector{FEInference.CPAmplitudeChannel}
  JET.@test_opt target_modules=(FloquetExpansions,) FEInference.reconstruct_cp_amplitude_channels(
    [transported]
  )

  generator = @inferred FEInference.reconstruct_cp_effective_generator(coherent, rows)
  @test generator isa Liouvillian
  JET.@test_opt target_modules=(FloquetExpansions,) FEInference.reconstruct_cp_effective_generator(
    coherent, rows
  )
end

using Test
using FloquetExpansions

include(joinpath(@__DIR__, "helpers", "magnus_log_recurrence.jl"))

const FE_MLRT = FloquetExpansions

@testset "Magnus logarithm recurrence matches exact canonical log through order six" begin
  R = Rational{Int}
  support = [-1, 0, 1]
  components = Dict(-1 => R[0 1; 2 -1], 0 => R[1 2; -1 0], 1 => R[2 -1; 1 1])
  zero_component = zeros(R, 2, 2)
  order = 6

  plan = FE_MLRT.compile_bloch_projection_plan(support, order)
  bloch = FE_MLRT.evaluate_bloch_projection_plan(
    plan, components; product=(*), inverse_weight=harmonic -> 1 // harmonic, zero_component
  )
  direct = FE_MLRT.bloch_van_vleck_reconstruction(
    plan, bloch; product=(*), zero_component
  )
  connected, counts = connected_log_magnus_recurrence(
    direct.normalized_embedding; product=(*), zero_component
  )

  @test connected == direct.log_embedding
  @test counts.harmonic_products > 0
end

@test mlr_bernoulli_numbers(6) == Rational{Int}[1, -1 // 2, 1 // 6, 0, -1 // 30, 0, 1 // 42]

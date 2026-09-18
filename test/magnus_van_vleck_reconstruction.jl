using Test
using FloquetExpansions

include(joinpath(@__DIR__, "helpers", "magnus_van_vleck_reconstruction.jl"))

const FE_MVVRT = FloquetExpansions

@testset "direct Magnus Van Vleck reconstruction matches exact reference through order six" begin
  R = Rational{Int}
  support = [-1, 0, 1]
  components = Dict(-1 => R[0 1; 2 -1], 0 => R[1 2; -1 0], 1 => R[2 -1; 1 1])
  zero_component = zeros(R, 2, 2)
  order = 6

  plan = FE_MVVRT.compile_bloch_projection_plan(support, order)
  bloch = FE_MVVRT.evaluate_bloch_projection_plan(
    plan, components; product=(*), inverse_weight=harmonic -> 1 // harmonic, zero_component
  )
  direct = FE_MVVRT.bloch_van_vleck_reconstruction(
    plan, bloch; product=(*), zero_component
  )
  magnus = magnus_van_vleck_reconstruction(plan, bloch; product=(*), zero_component)

  @test magnus.static_factor == direct.static_factor
  @test magnus.inverse_static_factor == direct.inverse_static_factor
  @test magnus.normalized_embedding == direct.normalized_embedding
  @test magnus.log_embedding == direct.log_embedding
  @test magnus.effective == direct.effective
  @test all(
    iszero(get(magnus.log_embedding[n], plan.zero_harmonic, zero_component)) for
    n in eachindex(magnus.log_embedding)
  )
end

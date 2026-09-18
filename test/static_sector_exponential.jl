using Test
using FloquetExpansions

include(joinpath(@__DIR__, "helpers", "connected_van_vleck_reconstruction.jl"))
include(joinpath(@__DIR__, "helpers", "static_sector_exponential.jl"))

const FE_SSE = FloquetExpansions

@testset "static-sector exponential reproduces full periodic convolution" begin
  R = Rational{Int}
  support = [-1, 0, 1]
  components = Dict(-1 => R[0 1; 2 -1], 0 => R[1 2; -1 0], 1 => R[2 -1; 1 1])
  zero_component = zeros(R, 2, 2)
  order = 6

  projection_plan = FE_SSE.compile_bloch_projection_plan(support, order)
  bloch = FE_SSE.evaluate_bloch_projection_plan(
    projection_plan,
    components;
    product=(*),
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component,
  )
  log_plan = compile_lyndon_log_evaluation_plan(support, order)
  log_embedding = evaluate_lyndon_log_plan(
    log_plan, components; product=(*), zero_component
  )
  identity_component = one(first(bloch.effective))

  full, full_counts = connected_static_factor(
    log_embedding,
    projection_plan.zero_harmonic,
    identity_component,
    zero_component;
    product=(*),
  )
  static_plan = compile_static_sector_exp_plan(log_plan, projection_plan.zero_harmonic)
  compiled = evaluate_static_sector_exp_plan(
    static_plan, log_embedding, identity_component, zero_component; product=(*)
  )

  @test compiled == full
  @test static_plan.product_count <= full_counts.harmonic_products
  @info "static-sector exponential profile" order static_plan.product_count full_products=full_counts.harmonic_products nodes=length(
    static_plan.nodes
  )
end

using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA_BVV = SecondQuantizedAlgebra

include(joinpath(@__DIR__, "helpers", "bloch_van_vleck_projection.jl"))

hamiltonian_product_projection(left::SQA_BVV.QAdd, right::SQA_BVV.QAdd) = left * right
liouvillian_product_projection(left::Liouvillian, right::Liouvillian) = compose(left, right)

projection_vanishes(value) = iszero(SQA_BVV.simplify(value))
function projection_vanishes(generator::PeriodicGenerator)
  return all(projection_vanishes(generator[harmonic]) for harmonic in keys(generator))
end

space_projection = NLevelSpace(:bloch_van_vleck_projection, 3)
σ11_projection = Transition(space_projection, :σ, 1, 1)
σ22_projection = Transition(space_projection, :σ, 2, 2)
σ33_projection = Transition(space_projection, :σ, 3, 3)
σ12_projection = Transition(space_projection, :σ, 1, 2)
σ23_projection = Transition(space_projection, :σ, 2, 3)
σ31_projection = Transition(space_projection, :σ, 3, 1)
σ13_projection = Transition(space_projection, :σ, 1, 3)
σ21_projection = Transition(space_projection, :σ, 2, 1)
σ32_projection = Transition(space_projection, :σ, 3, 2)
@variables ω_projection::Real

H0_projection =
  2 * σ11_projection - σ22_projection + (3 // 2) * σ33_projection +
  σ12_projection + σ12_projection'
H1_projection = σ12_projection + 2 * σ23_projection + im * σ31_projection
H2_projection = 2 * σ13_projection - σ21_projection + im * σ32_projection
H_projection = PeriodicGenerator(
  Dict(
    0 => H0_projection,
    1 => H1_projection,
    -1 => H1_projection',
    2 => H2_projection,
    -2 => H2_projection',
  ),
  ω_projection,
)

@testset "production Bloch core reconstructs canonical Hamiltonian Van Vleck" begin
  order = 3
  bloch = evaluate_periodic_bloch_projection(
    H_projection,
    order;
    product=hamiltonian_product_projection,
    inverse_weight=harmonic -> 1 // harmonic,
    simplifier=SQA_BVV.simplify,
  )
  converted = bloch_van_vleck_projection(
    bloch, H_projection; product=hamiltonian_product_projection, simplifier=SQA_BVV.simplify
  )
  vv = floquet_expansion(H_projection, VanVleck(), order)

  @test converted.counts.log_products == expected_projection_mercator_products(order - 1)
  @test all(projection_vanishes(time_average(term)) for term in converted.log_embedding)
  for n in 1:(order - 1)
    @test projection_vanishes(
      converted.effective[n + 1] - ω_projection^n * effective_component(vv, n)
    )
    @test projection_vanishes(
      im * converted.log_embedding[n] - ω_projection^n * micromotion(vv, n)
    )
  end

  replayed_log, log_counts = replay_mercator_log(
    converted.normalized_embedding,
    H_projection;
    product=hamiltonian_product_projection,
    simplifier=SQA_BVV.simplify,
  )
  @test log_counts.log_products == converted.counts.log_products
  @test all(
    projection_vanishes(replayed_log[n] - converted.log_embedding[n]) for
    n in eachindex(replayed_log)
  )

  replayed_effective, _ = replay_static_similarity(
    bloch.effective,
    converted.static_factor;
    product=hamiltonian_product_projection,
    simplifier=SQA_BVV.simplify,
  )
  @test all(
    projection_vanishes(replayed_effective[n] - converted.effective[n]) for
    n in eachindex(replayed_effective)
  )
end

L_projection = PeriodicGenerator(
  Dict(
    0 =>
      hamiltonian_action(σ11_projection - σ33_projection) +
      dissipator(σ12_projection + σ23_projection),
    1 => hamiltonian_action(σ12_projection + σ21_projection),
    -1 => dissipator(σ23_projection + σ31_projection),
    2 => hamiltonian_action(σ13_projection + σ31_projection),
    -2 => dissipator(σ12_projection - im * σ23_projection),
  ),
  ω_projection,
)

@testset "production Bloch core reconstructs canonical Liouvillian Van Vleck" begin
  order = 3
  bloch = evaluate_periodic_bloch_projection(
    L_projection,
    order;
    product=liouvillian_product_projection,
    inverse_weight=harmonic -> im // harmonic,
    simplifier=SQA_BVV.simplify,
  )
  converted = bloch_van_vleck_projection(
    bloch, L_projection; product=liouvillian_product_projection, simplifier=SQA_BVV.simplify
  )
  vv = floquet_expansion(L_projection, VanVleck(), order)

  @test converted.counts.log_products == expected_projection_mercator_products(order - 1)
  @test all(projection_vanishes(time_average(term)) for term in converted.log_embedding)
  for n in 1:(order - 1)
    @test projection_vanishes(
      converted.effective[n + 1] - ω_projection^n * effective_component(vv, n)
    )
    @test projection_vanishes(
      converted.log_embedding[n] - ω_projection^n * micromotion(vv, n)
    )
  end
end

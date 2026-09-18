using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const FE_BCVV = FloquetExpansions
const SQA_BCVV = SecondQuantizedAlgebra

function bcvv_embedding_equal(left, right, zero_component; simplifier=identity)
  harmonics = union(keys(left), keys(right))
  return all(
    iszero(
      simplifier(get(left, harmonic, zero_component) - get(right, harmonic, zero_component))
    ) for harmonic in harmonics
  )
end

@testset "connected Van Vleck matches dense exact reference through order six" begin
  R = Rational{Int}
  support = [-1, 0, 1]
  components = Dict(-1 => R[0 1; 2 -1], 0 => R[1 2; -1 0], 1 => R[2 -1; 1 1])
  zero_component = zeros(R, 2, 2)
  order = 6

  projection_plan = FE_BCVV.compile_bloch_projection_plan(support, order)
  bloch = FE_BCVV.evaluate_bloch_projection_plan(
    projection_plan,
    components;
    product=(*),
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component,
  )
  reference = FE_BCVV.bloch_van_vleck_reconstruction(
    projection_plan, bloch; product=(*), zero_component
  )
  connected_plan = FE_BCVV.compile_bloch_connected_van_vleck_plan(projection_plan)
  connected = FE_BCVV.evaluate_bloch_connected_van_vleck(
    projection_plan, bloch, components, connected_plan; product=(*), zero_component
  )

  @test connected.static_factor == reference.static_factor
  @test connected.inverse_static_factor == reference.inverse_static_factor
  @test connected.effective == reference.effective
  @test length(connected.log_embedding) == length(reference.log_embedding)
  @test all(
    bcvv_embedding_equal(
      connected.log_embedding[n], reference.log_embedding[n], zero_component
    ) for n in eachindex(reference.log_embedding)
  )
  @test all(
    !haskey(connected.log_embedding[n], projection_plan.zero_harmonic) for
    n in eachindex(connected.log_embedding)
  )
  @test length(connected_plan.log.brackets) == 66
  @test FE_BCVV.bloch_connected_log_products(connected_plan.log) == 132
  @test connected_plan.static.product_count == 52
end

@testset "connected Van Vleck matches SQA reference" begin
  space = PauliSpace(:bloch_connected_vv)
  sigma_x = Pauli(space, :sigma, 1)
  sigma_y = Pauli(space, :sigma, 2)
  sigma_z = Pauli(space, :sigma, 3)
  @variables omega_bcvv::Real

  drive = sigma_x + im * sigma_y
  generator = PeriodicGenerator(
    Dict(0 => 1 * sigma_z, 1 => drive, -1 => drive'), omega_bcvv
  )
  components = getfield(generator, :components)
  zero_component = getfield(generator, :zero_component)
  product(left, right) = left * right
  order = 4

  projection_plan = FE_BCVV.compile_bloch_projection_plan(collect(keys(generator)), order)
  bloch = FE_BCVV.evaluate_bloch_projection_plan(
    projection_plan,
    components;
    product,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component,
    simplifier=SQA_BCVV.simplify,
  )
  reference = FE_BCVV.bloch_van_vleck_reconstruction(
    projection_plan, bloch; product, zero_component, simplifier=SQA_BCVV.simplify
  )
  connected_plan = FE_BCVV.compile_bloch_connected_van_vleck_plan(projection_plan)
  connected = FE_BCVV.evaluate_bloch_connected_van_vleck(
    projection_plan,
    bloch,
    components,
    connected_plan;
    product,
    zero_component,
    simplifier=SQA_BCVV.simplify,
  )

  @test all(
    iszero(SQA_BCVV.simplify(connected.static_factor[n] - reference.static_factor[n])) for
    n in eachindex(reference.static_factor)
  )
  @test all(
    iszero(SQA_BCVV.simplify(connected.effective[n] - reference.effective[n])) for
    n in eachindex(reference.effective)
  )
  @test all(
    bcvv_embedding_equal(
      connected.log_embedding[n],
      reference.log_embedding[n],
      zero_component;
      simplifier=SQA_BCVV.simplify,
    ) for n in eachindex(reference.log_embedding)
  )
end

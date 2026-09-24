using Test
using LinearAlgebra: norm
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const FE_BVVI = FloquetExpansions
const SQA_BVVI = SecondQuantizedAlgebra

struct BVVIHarmonic2D
  n1::Int
  n2::Int
end

function Base.:+(left::BVVIHarmonic2D, right::BVVIHarmonic2D)
  return BVVIHarmonic2D(left.n1 + right.n1, left.n2 + right.n2)
end
Base.zero(::BVVIHarmonic2D) = BVVIHarmonic2D(0, 0)
Base.iszero(harmonic::BVVIHarmonic2D) = iszero(harmonic.n1) && iszero(harmonic.n2)

function bvvi_vanishes(value)
  return iszero(SQA_BVVI.simplify(value))
end

function bvvi_vanishes(generator::PeriodicGenerator)
  return all(bvvi_vanishes(generator[harmonic]) for harmonic in keys(generator))
end

function bvvi_periodic(embedding, template)
  return PeriodicGenerator(
    Dict(harmonic => value for (harmonic, value) in embedding),
    template.wd,
    template.zero_component,
  )
end

function bvvi_projection(generator, order; product, inverse_weight, simplifier=identity)
  components = getfield(generator, :components)
  plan = FE_BVVI.compile_bloch_projection_plan(collect(keys(generator)), order)
  result = FE_BVVI.evaluate_bloch_projection_plan(
    plan,
    components;
    product,
    inverse_weight,
    zero_component=getfield(generator, :zero_component),
    simplifier,
  )
  return plan, result
end

@testset "internal Bloch reconstruction matches Hamiltonian Van Vleck" begin
  space = PauliSpace(:bloch_vv_internal_hamiltonian)
  σx = Pauli(space, :σ, 1)
  σy = Pauli(space, :σ, 2)
  σz = Pauli(space, :σ, 3)
  @variables ω_bvvi::Real

  H1 = σx + im * σy
  H = PeriodicGenerator(Dict(0 => 1 * σz, 1 => H1, -1 => H1'), ω_bvvi)
  product(left, right) = left * right
  order = 3

  plan, bloch = bvvi_projection(
    H,
    order;
    product,
    inverse_weight=harmonic -> 1 // harmonic,
    simplifier=SQA_BVVI.simplify,
  )
  converted = FE_BVVI.bloch_van_vleck_reconstruction(
    plan,
    bloch;
    product,
    zero_component=getfield(H, :zero_component),
    simplifier=SQA_BVVI.simplify,
  )
  vv = floquet_expansion(H, VanVleck(), order)

  @test converted.counts.log_products ==
    FE_BVVI.bloch_van_vleck_mercator_products(order - 1)
  @test all(!haskey(term, plan.zero_harmonic) for term in converted.log_embedding)
  for n in 1:(order - 1)
    @test bvvi_vanishes(converted.effective[n + 1] - ω_bvvi^n * effective_component(vv, n))
    kick = im * bvvi_periodic(converted.log_embedding[n], H)
    @test bvvi_vanishes(kick - ω_bvvi^n * micromotion(vv, n))
  end
end

@testset "internal Bloch reconstruction matches generic Liouvillian Van Vleck" begin
  space = PauliSpace(:bloch_vv_internal_liouvillian)
  σx = Pauli(space, :σ, 1)
  σy = Pauli(space, :σ, 2)
  σz = Pauli(space, :σ, 3)
  @variables ω_bvvi_map::Real

  L = PeriodicGenerator(
    Dict(
      0 => hamiltonian_action(σz) + dissipator(σx),
      1 => hamiltonian_action(σx + σy),
      -1 => dissipator(σy + σz),
    ),
    ω_bvvi_map,
  )
  product(left, right) = compose(left, right)
  order = 3

  plan, bloch = bvvi_projection(
    L,
    order;
    product,
    inverse_weight=harmonic -> im // harmonic,
    simplifier=SQA_BVVI.simplify,
  )
  converted = FE_BVVI.bloch_van_vleck_reconstruction(
    plan,
    bloch;
    product,
    zero_component=getfield(L, :zero_component),
    simplifier=SQA_BVVI.simplify,
  )
  vv = floquet_expansion(
    L, VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))), order
  )

  @test all(!haskey(term, plan.zero_harmonic) for term in converted.log_embedding)
  for n in 1:(order - 1)
    @test bvvi_vanishes(
      converted.effective[n + 1] - ω_bvvi_map^n * effective_component(vv, n)
    )
    kick = bvvi_periodic(converted.log_embedding[n], L)
    @test bvvi_vanishes(kick - ω_bvvi_map^n * micromotion(vv, n))
  end
end

@testset "internal Bloch reconstruction preserves generic additive harmonics" begin
  zero_harmonic = BVVIHarmonic2D(0, 0)
  h10 = BVVIHarmonic2D(1, 0)
  hm10 = BVVIHarmonic2D(-1, 0)
  h01 = BVVIHarmonic2D(0, 1)
  h0m1 = BVVIHarmonic2D(0, -1)
  support = [zero_harmonic, h10, hm10, h01, h0m1]
  components = Dict{BVVIHarmonic2D,Matrix{ComplexF64}}(
    zero_harmonic => ComplexF64[0.2 0.3im; -0.1 0.4],
    h10 => ComplexF64[0.1 0.7; 0.2im -0.3],
    hm10 => ComplexF64[0.3im -0.2; 0.5 0.1],
    h01 => ComplexF64[-0.2 0.4im; 0.6 0.3],
    h0m1 => ComplexF64[0.1 -0.5; 0.2 0.4im],
  )
  product(left, right) = left * right
  inverse_weight(harmonic) = inv(harmonic.n1 + sqrt(2) * harmonic.n2)
  order = 4
  zero_component = zeros(ComplexF64, 2, 2)

  plan = FE_BVVI.compile_bloch_projection_plan(support, order, zero_harmonic)
  bloch = FE_BVVI.evaluate_bloch_projection_plan(
    plan, components; product, inverse_weight, zero_component
  )
  converted = FE_BVVI.bloch_van_vleck_reconstruction(plan, bloch; product, zero_component)

  @test length(converted.effective) == order
  @test length(converted.log_embedding) == order - 1
  @test converted.counts.log_products ==
    FE_BVVI.bloch_van_vleck_mercator_products(order - 1)
  @test all(!haskey(term, zero_harmonic) for term in converted.log_embedding)
  @test all(
    all(harmonic isa BVVIHarmonic2D for harmonic in keys(term)) for
    term in converted.log_embedding
  )
  @test all(isfinite, (norm(term) for term in converted.effective))
end

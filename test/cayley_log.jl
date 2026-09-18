using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA_CAYLEY = SecondQuantizedAlgebra

include(joinpath(@__DIR__, "helpers", "cayley_log.jl"))

function cayley_exact_equal(left::Dict{H,T}, right::Dict{H,T}) where {H,T}
  keys_union = union(keys(left), keys(right))
  return all(get(left, key, zero(first(values(left)))) == get(right, key, zero(first(values(right)))) for key in keys_union)
end

@testset "Cayley log matches Mercator in a noncommutative matrix algebra" begin
  R = Rational{Int}
  A = Matrix{R}[
    [0 1; 2 0],
    [1 1; 0 -1],
    [2 -1; 1 1],
    [0 2; -1 1],
    [1 0; 3 -2],
    [2 1; -2 0],
  ]
  B = Matrix{R}[
    [1 0; 1 -1],
    [0 1; -1 2],
    [2 1; 0 -1],
    [1 -2; 1 0],
    [0 1; 2 1],
    [1 2; -1 1],
  ]
  normalized = [Dict(-1 => A[n], 1 => B[n]) for n in 1:6]
  product(left, right) = left * right

  mercator, mercator_counts = mercator_log_series(normalized; product)
  cayley, cayley_counts = cayley_log_series(normalized; product)

  @test cayley == mercator
  @test mercator_counts.power_products == 35
  @test cayley_counts.cayley_products == 15
  @test cayley_counts.power_products == 28
  @test cayley_counts.cayley_products + cayley_counts.power_products == 43
end

@testset "Cayley log matches the certified SQA canonical logarithm" begin
  space = PauliSpace(:cayley_log_pauli)
  σx = Pauli(space, :σ, 1)
  σy = Pauli(space, :σ, 2)
  σz = Pauli(space, :σ, 3)
  @variables ω_cayley::Real

  H1 = σx + im * σy
  H = PeriodicGenerator(Dict(0 => 1 * σz, 1 => H1, -1 => H1'), ω_cayley)
  components = getfield(H, :components)
  order = 5
  plan = FloquetExpansions.compile_bloch_projection_plan(collect(keys(H)), order)
  bloch = FloquetExpansions.evaluate_bloch_projection_plan(
    plan,
    components;
    product=*,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component=getfield(H, :zero_component),
    simplifier=SQA_CAYLEY.simplify,
  )
  converted = FloquetExpansions.bloch_van_vleck_reconstruction(
    plan,
    bloch;
    product=*,
    zero_component=getfield(H, :zero_component),
    simplifier=SQA_CAYLEY.simplify,
  )
  cayley, _ = cayley_log_series(
    converted.normalized_embedding; product=*, simplifier=SQA_CAYLEY.simplify
  )

  @test length(cayley) == length(converted.log_embedding)
  for n in eachindex(cayley)
    harmonics = union(keys(cayley[n]), keys(converted.log_embedding[n]))
    for harmonic in harmonics
      left = get(cayley[n], harmonic, getfield(H, :zero_component))
      right = get(converted.log_embedding[n], harmonic, getfield(H, :zero_component))
      @test iszero(SQA_CAYLEY.simplify(left - right))
    end
  end
end

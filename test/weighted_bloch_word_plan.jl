using Test
using FloquetExpansions
using LinearAlgebra: norm
using Symbolics: @variables

include(joinpath(@__DIR__, "helpers", "bloch_feshbach_reference.jl"))
include(joinpath(@__DIR__, "helpers", "harmonic_word_plan.jl"))
include(joinpath(@__DIR__, "helpers", "weighted_bloch_word_plan.jl"))

flow_word_set(words) = Set{Tuple}(word.harmonics for word in words)

function evaluate_flow_word(generator, flow_word, product)
  isempty(flow_word) && throw(ArgumentError("flow word must not be empty"))
  value = generator[flow_word[end]]
  for index in (length(flow_word) - 1):-1:1
    value = product(value, generator[flow_word[index]])
  end
  return value
end

function evaluate_word_sum(terms, generator, product, zero_harmonic)
  value = zero(generator[zero_harmonic])
  for (word, coefficient) in terms
    value += coefficient * evaluate_flow_word(generator, word, product)
  end
  return value
end

function evaluate_wave_component(
  plan, generator, product, order_index, harmonic, zero_harmonic
)
  terms = get(plan.wave[order_index], harmonic, nothing)
  isnothing(terms) && return zero(generator[zero_harmonic])
  return evaluate_word_sum(terms, generator, product, zero_harmonic)
end

@testset "weighted Bloch plan attaches exact primitive and fold coefficients" begin
  support = [-2, -1, 0, 1, 2]
  inverse_weight = harmonic -> 1 // harmonic
  plan = compile_bloch_word_plan(support, 4; inverse_weight)

  for degree in 1:4
    topology = compile_harmonic_word_plan(support, degree)
    weighted_words = Set(keys(plan.effective[degree]))
    @test issubset(weighted_words, flow_word_set(topology.closed))
  end

  for degree in 2:4
    topology = compile_harmonic_word_plan(support, degree)
    for word in first_return_words(topology)
      @test haskey(plan.effective[degree], word.harmonics)
      @test plan.effective[degree][word.harmonics] ==
        primitive_bloch_word_weight(word.harmonics, inverse_weight, 0)
    end
  end

  for harmonic in filter(!iszero, support)
    folded_word = (0, harmonic, -harmonic)
    @test plan.effective[3][folded_word] == -(inverse_weight(harmonic)^2)
  end
  @test plan.effective[3][(0, 2, -2)] != inverse_weight(2)

  @test plan.counts.generator_product_terms > 0
  @test plan.counts.folded_counterterms > 0
end

@testset "weighted Bloch plan preserves Hamiltonian and Liouvillian homological phases" begin
  support = [-2, -1, 0, 1, 2]
  hamiltonian_weight = harmonic -> 1 // harmonic
  liouvillian_weight = harmonic -> im // harmonic
  hamiltonian_plan = compile_bloch_word_plan(support, 4; inverse_weight=hamiltonian_weight)
  liouvillian_plan = compile_bloch_word_plan(support, 4; inverse_weight=liouvillian_weight)

  for degree in 2:4
    topology = compile_harmonic_word_plan(support, degree)
    for word in first_return_words(topology)
      @test hamiltonian_plan.effective[degree][word.harmonics] ==
        primitive_bloch_word_weight(word.harmonics, hamiltonian_weight, 0)
      @test liouvillian_plan.effective[degree][word.harmonics] ==
        primitive_bloch_word_weight(word.harmonics, liouvillian_weight, 0)
    end
  end
end

struct WeightedTestHarmonicIndex
  n1::Int
  n2::Int
end

Base.zero(::Type{WeightedTestHarmonicIndex}) = WeightedTestHarmonicIndex(0, 0)
Base.zero(::WeightedTestHarmonicIndex) = zero(WeightedTestHarmonicIndex)
function Base.iszero(index::WeightedTestHarmonicIndex)
  return iszero(index.n1) && iszero(index.n2)
end
function Base.:+(left::WeightedTestHarmonicIndex, right::WeightedTestHarmonicIndex)
  return WeightedTestHarmonicIndex(left.n1 + right.n1, left.n2 + right.n2)
end
Base.:-(index::WeightedTestHarmonicIndex) = WeightedTestHarmonicIndex(-index.n1, -index.n2)

@testset "weighted Bloch plan remains generic over multidimensional harmonic labels" begin
  zero_index = zero(WeightedTestHarmonicIndex)
  e1 = WeightedTestHarmonicIndex(1, 0)
  e2 = WeightedTestHarmonicIndex(0, 1)
  support = [zero_index, e1, -e1, e2, -e2]
  inverse_weight = index -> 1 // (7 * index.n1 + 11 * index.n2)
  plan = compile_bloch_word_plan(support, 4; inverse_weight, zero_harmonic=zero_index)

  for degree in 2:4
    topology = compile_harmonic_word_plan(support, degree; zero_harmonic=zero_index)
    for word in first_return_words(topology)
      @test haskey(plan.effective[degree], word.harmonics)
      @test plan.effective[degree][word.harmonics] ==
        primitive_bloch_word_weight(word.harmonics, inverse_weight, zero_index)
    end
  end
end

@testset "weighted Bloch plan reproduces the certified recurrence before operator algebra" begin
  @variables w::Real
  matrix_generator = PeriodicGenerator(
    Dict(
      0 => ComplexF64[0.4 1.0 + 0.2im; -0.3 + 0.1im -0.7],
      1 => ComplexF64[0.2 + 0.1im -0.6; 0.5im 0.8 - 0.2im],
      -1 => ComplexF64[-0.1 + 0.4im 0.3; 0.7 - 0.2im 0.5im],
      2 => ComplexF64[0.6 -0.2im; 0.9 + 0.1im -0.4],
      -2 => ComplexF64[-0.5im 0.8; -0.1 + 0.3im 0.2],
    ),
    w,
  )
  product = (left, right) -> left * right
  zero_harmonic = 0
  support = collect(keys(matrix_generator))
  order = 6

  for inverse_weight in (harmonic -> 1.0 / harmonic, harmonic -> im / harmonic)
    reference = bloch_reference(matrix_generator, order; product, inverse_weight)
    plan = compile_bloch_word_plan(support, order; inverse_weight, zero_harmonic)

    @test length(plan.effective) == length(reference.effective)
    @test length(plan.wave) == length(reference.wave)

    for n in eachindex(reference.effective)
      planned = evaluate_word_sum(
        plan.effective[n], matrix_generator, product, zero_harmonic
      )
      @test norm(planned - reference.effective[n]) <= 5.0e-11
      @test all(length(word) == n for word in keys(plan.effective[n]))
      @test all(iszero(sum(word)) for word in keys(plan.effective[n]))
    end

    for n in eachindex(reference.wave)
      harmonics = union(keys(reference.wave[n]), keys(plan.wave[n]))
      for harmonic in harmonics
        planned = evaluate_wave_component(
          plan, matrix_generator, product, n, harmonic, zero_harmonic
        )
        @test norm(planned - reference.wave[n][harmonic]) <= 5.0e-11
      end
      @test all(
        length(word) == n && sum(word) == harmonic for (harmonic, terms) in plan.wave[n] for
        word in keys(terms)
      )
    end
  end
end

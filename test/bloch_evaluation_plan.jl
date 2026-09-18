using Test
using FloquetExpansions
using LinearAlgebra: norm
using Symbolics: @variables

include(joinpath(@__DIR__, "helpers", "bloch_feshbach_reference.jl"))
include(joinpath(@__DIR__, "helpers", "harmonic_word_plan.jl"))
include(joinpath(@__DIR__, "helpers", "weighted_bloch_word_plan.jl"))
include(joinpath(@__DIR__, "helpers", "bloch_evaluation_plan.jl"))

function matrix_components(generator)
  return Dict(harmonic => generator[harmonic] for harmonic in keys(generator))
end

function weighted_word_term_count(plan)
  effective_terms = sum(length, plan.effective)
  wave_terms = sum(
    length(terms) for harmonic_map in plan.wave for terms in values(harmonic_map)
  )
  return effective_terms + wave_terms
end

@testset "sparse Bloch evaluation plan reproduces the certified dense recurrence" begin
  @variables w::Real
  generator = PeriodicGenerator(
    Dict(
      0 => ComplexF64[0.4 1.0 + 0.2im; -0.3 + 0.1im -0.7],
      1 => ComplexF64[0.2 + 0.1im -0.6; 0.5im 0.8 - 0.2im],
      -1 => ComplexF64[-0.1 + 0.4im 0.3; 0.7 - 0.2im 0.5im],
      2 => ComplexF64[0.6 -0.2im; 0.9 + 0.1im -0.4],
      -2 => ComplexF64[-0.5im 0.8; -0.1 + 0.3im 0.2],
    ),
    w,
  )
  components = matrix_components(generator)
  support = collect(keys(generator))
  zero_component = zero(generator[0])
  product = (left, right) -> left * right
  order = 6
  plan = compile_bloch_evaluation_plan(support, order)

  for inverse_weight in (harmonic -> 1.0 / harmonic, harmonic -> im / harmonic)
    reference = bloch_reference(generator, order; product, inverse_weight)
    evaluated = evaluate_bloch_evaluation_plan(
      plan, components; product, inverse_weight, zero_component
    )

    @test length(evaluated.effective) == length(reference.effective)
    @test length(evaluated.wave) == length(reference.wave)
    for n in eachindex(reference.effective)
      @test norm(evaluated.effective[n] - reference.effective[n]) <= 5.0e-11
    end
    for n in eachindex(reference.wave)
      @test Set(keys(evaluated.wave[n])) == Set(keys(reference.wave[n]))
      for harmonic in keys(reference.wave[n])
        @test norm(evaluated.wave[n][harmonic] - reference.wave[n][harmonic]) <= 5.0e-11
      end
    end

    @test plan.counts.generator_products + plan.counts.fold_products ==
      reference.counts.component_products
  end
end

@testset "sparse recurrence states compact the expanded weighted-word reference" begin
  support = [-2, -1, 0, 1, 2]
  order = 6
  evaluation_plan = compile_bloch_evaluation_plan(support, order)
  word_plan = compile_bloch_word_plan(
    support, order; inverse_weight=harmonic -> 1 // harmonic
  )

  recurrence_states =
    evaluation_plan.counts.residual_nodes +
    evaluation_plan.counts.wave_nodes +
    evaluation_plan.counts.effective_nodes
  expanded_words = weighted_word_term_count(word_plan)

  @test evaluation_plan.counts.residual_nodes == 85
  @test evaluation_plan.counts.wave_nodes == 60
  @test evaluation_plan.counts.effective_nodes == 6
  @test evaluation_plan.counts.generator_products == 300
  @test evaluation_plan.counts.fold_products == 140
  @test recurrence_states == 151
  @test recurrence_states < expanded_words
  @test all(
    !iszero(harmonic) for support_n in evaluation_plan.wave_support for
    harmonic in support_n
  )
end

struct EvaluationTestHarmonicIndex
  n1::Int
  n2::Int
end

Base.zero(::Type{EvaluationTestHarmonicIndex}) = EvaluationTestHarmonicIndex(0, 0)
Base.zero(::EvaluationTestHarmonicIndex) = zero(EvaluationTestHarmonicIndex)
function Base.iszero(index::EvaluationTestHarmonicIndex)
  return iszero(index.n1) && iszero(index.n2)
end
function Base.:+(left::EvaluationTestHarmonicIndex, right::EvaluationTestHarmonicIndex)
  return EvaluationTestHarmonicIndex(left.n1 + right.n1, left.n2 + right.n2)
end
function Base.:-(index::EvaluationTestHarmonicIndex)
  return EvaluationTestHarmonicIndex(-index.n1, -index.n2)
end

@testset "sparse evaluation-plan geometry is not scalar-frequency specific" begin
  zero_index = zero(EvaluationTestHarmonicIndex)
  e1 = EvaluationTestHarmonicIndex(1, 0)
  e2 = EvaluationTestHarmonicIndex(0, 1)
  support = [zero_index, e1, -e1, e2, -e2]
  plan = compile_bloch_evaluation_plan(support, 5; zero_harmonic=zero_index)

  @test plan.zero_harmonic == zero_index
  @test plan.effective_active[1]
  @test all(!iszero(harmonic) for support_n in plan.wave_support for harmonic in support_n)
  @test plan.counts.residual_nodes == sum(length, plan.residuals)
  @test plan.counts.wave_nodes == sum(length, plan.wave_support)
  @test plan.counts.generator_products > 0
  @test plan.counts.fold_products > 0
end

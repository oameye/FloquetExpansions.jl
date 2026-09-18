using Test
using FloquetExpansions

include(joinpath(@__DIR__, "helpers", "connected_word_log.jl"))

function cwl_compare_embeddings(left, right, zero_component)
  harmonics = union(keys(left), keys(right))
  return all(
    get(left, harmonic, zero_component) == get(right, harmonic, zero_component) for
    harmonic in harmonics
  )
end

@testset "formal harmonic-word reconstruction matches dense noncommuting algebra" begin
  R = Rational{Int}
  support = [-1, 0, 1]
  components = Dict(-1 => R[0 1; 2 -1], 0 => R[1 2; -1 0], 1 => R[2 -1; 1 1])
  zero_component = zeros(R, 2, 2)
  order = 6
  plan = FloquetExpansions.compile_bloch_projection_plan(support, order)
  direct_bloch = FloquetExpansions.evaluate_bloch_projection_plan(
    plan, components; product=(*), inverse_weight=harmonic -> 1 // harmonic, zero_component
  )
  direct = FloquetExpansions.bloch_van_vleck_reconstruction(
    plan, direct_bloch; product=(*), zero_component
  )

  word_plan, word_bloch, word = connected_word_reconstruction(support, order)
  @test word_plan.input_support == plan.input_support
  @test length(word_bloch.effective) == length(direct_bloch.effective)

  for n in eachindex(word.effective)
    evaluated = evaluate_harmonic_word_polynomial(
      word.effective[n], components; product=(*), zero_component
    )
    @test evaluated == direct.effective[n]
  end

  for n in eachindex(word.log_embedding)
    evaluated = evaluate_harmonic_word_embedding(
      word.log_embedding[n], components; product=(*), zero_component
    )
    @test cwl_compare_embeddings(evaluated, direct.log_embedding[n], zero_component)
  end
end

@testset "connected-word structural profile" begin
  support = [-1, 0, 1]
  for order in 2:6
    _, bloch, converted = connected_word_reconstruction(support, order)
    wave_words = harmonic_word_count(bloch.wave)
    normalized_words = harmonic_word_count(converted.normalized_embedding)
    log_words = harmonic_word_count(converted.log_embedding)
    prefix_products = harmonic_word_prefix_products(converted.log_embedding)
    println(
      "CONNECTED_WORD_PROFILE ",
      "order=$(order) ",
      "wave_words=$(wave_words) ",
      "normalized_words=$(normalized_words) ",
      "log_words=$(log_words) ",
      "prefix_products=$(prefix_products) ",
      "reconstruction_harmonic_products=$(converted.counts.harmonic_products)",
    )
    @test log_words > 0
    @test prefix_products >= 0
  end
end

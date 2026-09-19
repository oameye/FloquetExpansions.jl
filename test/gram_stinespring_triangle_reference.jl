using Test
using LinearAlgebra: I, kron

include(joinpath(@__DIR__, "helpers", "gram_stinespring_triangle_reference.jl"))

const GramExact = Complex{Rational{Int}}

function gram_exact_identity(dimension::Int)
  return Matrix{GramExact}(I, dimension, dimension)
end

function gram_recycling_super(left::Matrix{GramExact}, right::Matrix{GramExact})
  return kron(conj.(right), left)
end

function gram_all_period_coefficients_zero(polynomial::FormalPeriodMatrix)
  return all(iszero, polynomial.coefficients)
end

@testset "physical amplitude triangle is exactly a finite Q/R word-length cutoff" begin
  identity_component = gram_exact_identity(2)
  jump_a = gram_qr_jump(:a, 1, GramExact[1 1; 0 1])
  jump_b = gram_qr_jump(:b, -1, GramExact[0 1; 1 0])
  drift = gram_qr_drift(0, GramExact[1 0; 1 -1])

  for generator_order in 0:3
    words = gram_stinespring_triangle_words(
      [jump_a, jump_b], [drift], generator_order, identity_component
    )
    @test length(words) == gram_stinespring_triangle_count(2, 1, generator_order)
    @test length(words) == sum(3^length for length in 0:(generator_order + 1))
    @test all(word.output_number + word.drift_order <= generator_order + 1 for word in words)
  end

  words = gram_stinespring_triangle_words([jump_a, jump_b], [drift], 2, identity_component)
  for q in 0:3, d in 0:(3 - q)
    observed = count(word -> word.output_number == q && word.drift_order == d, words)
    @test observed == gram_stinespring_qd_count(2, 1, q, d)
  end
  @test all(word.output_number + word.drift_order <= 3 for word in words)
end

@testset "Gram metric is Hermitian without orthogonalizing physical output kernels" begin
  identity_component = gram_exact_identity(2)
  zero_component = zero(identity_component)
  jump_a = gram_qr_jump(:a, 1, GramExact[0 1; 1 0])
  jump_b = gram_qr_jump(:a, -1, GramExact[1 1; 0 -1])
  drift = gram_qr_drift(2, GramExact[1 0; 1 1])

  words = gram_stinespring_triangle_words(
    [jump_a, jump_b], [drift], 1, identity_component
  )
  metric = gram_stinespring_metric_series(words, 1, zero_component)

  for coefficient in metric, matrix in coefficient.coefficients
    @test matrix == adjoint(matrix)
  end
  @test formal_period_coefficient(metric[1], 0, zero_component) == identity_component
end

@testset "static scalar GKSL parent satisfies the triangle defect theorem exactly" begin
  identity_component = GramExact[1;;]
  zero_component = zero(identity_component)
  jump = gram_qr_jump(:loss, 0, GramExact[1;;])
  drift = gram_qr_drift(0, GramExact[-1 // 2;;])

  for generator_order in 0:2
    words = gram_stinespring_triangle_words(
      [jump], [drift], generator_order, identity_component
    )
    @test length(words) == sum(2^length for length in 0:(generator_order + 1))

    metric = gram_stinespring_metric_series(words, generator_order, zero_component)
    @test formal_period_coefficient(metric[1], 0, zero_component) == identity_component

    for channel_order in 1:(generator_order + 1)
      @test gram_all_period_coefficients_zero(metric[channel_order + 1])
    end

    defect = metric[generator_order + 3]
    @test !gram_all_period_coefficients_zero(defect)
  end
end

@testset "Gram-form channel and metric agree for a scalar physical system" begin
  identity_component = GramExact[1;;]
  zero_component = zero(identity_component)
  zero_superoperator = GramExact[0;;]
  jump = gram_qr_jump(:loss, 0, GramExact[1;;])
  drift = gram_qr_drift(0, GramExact[-1 // 2;;])

  words = gram_stinespring_triangle_words([jump], [drift], 1, identity_component)
  metric = gram_stinespring_metric_series(words, 1, zero_component)
  channel = gram_stinespring_channel_series(
    words, 1, zero_superoperator; paste=gram_recycling_super
  )

  @test length(channel) == length(metric)
  for channel_order in eachindex(metric)
    @test length(channel[channel_order].coefficients) ==
      length(metric[channel_order].coefficients)
    for period_degree in eachindex(metric[channel_order].coefficients)
      @test channel[channel_order].coefficients[period_degree] ==
        metric[channel_order].coefficients[period_degree]
    end
  end
end

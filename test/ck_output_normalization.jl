using Test
using FloquetExpansions
using LinearAlgebra: I

const CKNormalizationExact = Complex{Rational{Int}}
const ck_normalization_im = CKNormalizationExact(0 // 1, 1 // 1)

function ck_normalization_identity()
  return Matrix{CKNormalizationExact}(I, 2, 2)
end

function ck_normalization_polynomial(terms::Pair{Int,Matrix{CKNormalizationExact}}...)
  zero_component = zeros(CKNormalizationExact, 2, 2)
  return FloquetExpansions.CKOutputPairingPolynomial(Dict(terms), zero_component)
end

function ck_normalization_series(coefficients...)
  return FloquetExpansions.CKOutputPairingSeries(collect(coefficients))
end

function ck_normalization_vacuum_key()
  return FloquetExpansions.CKOutputKernelKey(
    Int[],
    FloquetExpansions.CKOutputBlock[],
    FloquetExpansions.CKOutputConstraint{Int}[],
    FloquetExpansions.CKOutputConstraint{Int}[],
  )
end

function ck_normalization_jump_key(harmonic::Int)
  return FloquetExpansions.CKOutputKernelKey(
    [1],
    [FloquetExpansions.CKOutputBlock(0, 1, 0, 1)],
    [FloquetExpansions.CKOutputConstraint(1, 1, harmonic)],
    FloquetExpansions.CKOutputConstraint{Int}[],
  )
end

function ck_normalization_amplitude(period_power::Int, terms...)
  zero_component = zeros(CKNormalizationExact, 2, 2)
  coefficients = FloquetExpansions.CKOutputKernel{Int,Matrix{CKNormalizationExact}}[
    FloquetExpansions.CKOutputKernel(
      Dict{FloquetExpansions.CKOutputKernelKey{Int},Matrix{CKNormalizationExact}}(),
      zero_component,
    ) for _ in 0:period_power
  ]
  coefficients[period_power + 1] = FloquetExpansions.CKOutputKernel(
    Dict{FloquetExpansions.CKOutputKernelKey{Int},Matrix{CKNormalizationExact}}(terms),
    zero_component,
  )
  return FloquetExpansions.CKOutputPeriodPolynomial(coefficients)
end

ck_normalization_metric_pair(left, right) = adjoint(left) * right

@testset "noncommutative CK right-normalization solves S M S = I coefficientwise" begin
  identity_component = ck_normalization_identity()
  zero_component = zero(identity_component)
  im = CKNormalizationExact(0 // 1, 1 // 1)

  m1_zero = CKNormalizationExact[1 1; 1 0]
  m1_one = CKNormalizationExact[0 im; -im 1]
  m2_zero = CKNormalizationExact[2 -im; im -1]
  m2_two = CKNormalizationExact[1 0; 0 -1]

  @test m1_zero * m2_zero != m2_zero * m1_zero

  metric = ck_normalization_series(
    ck_normalization_polynomial(0 => identity_component),
    ck_normalization_polynomial(0 => m1_zero, 1 => m1_one),
    ck_normalization_polynomial(0 => m2_zero, 2 => m2_two),
  )
  normalization = FloquetExpansions.ck_output_metric_inverse_sqrt_series(
    metric, 4, identity_component
  )
  normalized = FloquetExpansions.ck_output_normalized_metric_series(
    metric, normalization, 4, zero_component
  )

  @test normalized.coefficients[1].terms == Dict(0 => identity_component)
  @test all(isempty(coefficient.terms) for coefficient in normalized.coefficients[2:end])

  expected_s1 = ck_normalization_polynomial(
    0 => (-1 // 2) * m1_zero, 1 => (-1 // 2) * m1_one
  )
  @test normalization.coefficients[2].terms == expected_s1.terms

  for coefficient in normalization.coefficients, matrix in values(coefficient.terms)
    @test matrix == adjoint(matrix)
  end
end

@testset "CK normalization includes higher period powers without diagonalization" begin
  identity_component = ck_normalization_identity()
  zero_component = zero(identity_component)
  m1 = CKNormalizationExact[1 2; 2 -1]
  metric = ck_normalization_series(
    ck_normalization_polynomial(0 => identity_component),
    ck_normalization_polynomial(3 => m1),
  )

  normalization = FloquetExpansions.ck_output_metric_inverse_sqrt_series(
    metric, 3, identity_component
  )
  normalized = FloquetExpansions.ck_output_normalized_metric_series(
    metric, normalization, 3, zero_component
  )

  @test normalization.coefficients[2].terms == Dict(3 => (-1 // 2) * m1)
  @test haskey(normalization.coefficients[3].terms, 6)
  @test haskey(normalization.coefficients[4].terms, 9)
  @test normalized.coefficients[1].terms == Dict(0 => identity_component)
  @test all(isempty(coefficient.terms) for coefficient in normalized.coefficients[2:end])
end

@testset "right-normalized physical amplitudes reproduce S M S exactly" begin
  identity_component = ck_normalization_identity()
  zero_component = zero(identity_component)
  vacuum = ck_normalization_vacuum_key()
  jump = ck_normalization_jump_key(2)
  drift_value = CKNormalizationExact[0 1; 0 0]
  jump_value = CKNormalizationExact[1 0; 1 1]

  amplitudes = [
    ck_normalization_amplitude(0, vacuum => identity_component),
    ck_normalization_amplitude(1, vacuum => drift_value, jump => jump_value),
  ]
  metric = FloquetExpansions.ck_output_metric_series(
    amplitudes, 4, ck_normalization_im, zero_component, ck_normalization_metric_pair
  )
  normalization = FloquetExpansions.ck_output_metric_inverse_sqrt_series(
    metric, 4, identity_component
  )
  normalized_amplitudes = FloquetExpansions.ck_output_right_normalize_series(
    amplitudes, normalization, 4, zero_component
  )
  paired_normalized = FloquetExpansions.ck_output_metric_series(
    normalized_amplitudes,
    4,
    ck_normalization_im,
    zero_component,
    ck_normalization_metric_pair,
  )
  direct_normalized = FloquetExpansions.ck_output_normalized_metric_series(
    metric, normalization, 4, zero_component
  )

  @test paired_normalized.coefficients == direct_normalized.coefficients
  @test paired_normalized.coefficients[1].terms == Dict(0 => identity_component)
  @test all(
    isempty(coefficient.terms) for coefficient in paired_normalized.coefficients[2:end]
  )
  @test all(
    key.output_channels == [1] for coefficient in normalized_amplitudes for
    kernel in coefficient.coefficients for
    key in keys(kernel.terms) if !isempty(key.output_channels)
  )
end

@testset "CK normalization requires a unit leading metric" begin
  identity_component = ck_normalization_identity()
  bad_metric = ck_normalization_series(
    ck_normalization_polynomial(0 => 2 * identity_component)
  )
  @test_throws ArgumentError FloquetExpansions.ck_output_metric_inverse_sqrt_series(
    bad_metric, 2, identity_component
  )
end

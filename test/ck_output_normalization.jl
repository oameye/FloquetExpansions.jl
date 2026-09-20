using Test
using FloquetExpansions
using LinearAlgebra: I

const CKNormalizationExact = Complex{Rational{Int}}

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

@testset "CK normalization requires a unit leading metric" begin
  identity_component = ck_normalization_identity()
  bad_metric = ck_normalization_series(
    ck_normalization_polynomial(0 => 2 * identity_component)
  )
  @test_throws ArgumentError FloquetExpansions.ck_output_metric_inverse_sqrt_series(
    bad_metric, 2, identity_component
  )
end

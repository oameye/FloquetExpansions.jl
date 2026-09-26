using Test
using FloquetExpansions
using LinearAlgebra: I, kron

const CKNormalizationExact = Complex{Rational{Int}}
const ck_normalization_im = CKNormalizationExact(0 // 1, 1 // 1)

function ck_normalization_identity()
  return Matrix{CKNormalizationExact}(I, 2, 2)
end

function ck_normalization_polynomial(terms::Pair{Int,Matrix{CKNormalizationExact}}...)
  zero_component = zeros(CKNormalizationExact, 2, 2)
  return FloquetExpansions.CKOutputPairingPolynomial(Dict(terms), zero_component)
end

function ck_normalization_graded_polynomial(terms...)
  zero_component = zeros(CKNormalizationExact, 2, 2)
  graded = Dict{FloquetExpansions.CKOutputPairingGrade{Int},Matrix{CKNormalizationExact}}()
  for ((phase_harmonic, period_power), value) in terms
    graded[FloquetExpansions.CKOutputPairingGrade(phase_harmonic, period_power)] = value
  end
  return FloquetExpansions.CKOutputPairingPolynomial(graded, zero_component)
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

function ck_normalization_scalar_resolvent_key(harmonic::Int)
  return FloquetExpansions.CKOutputKernelKey(
    Int[],
    [FloquetExpansions.CKOutputBlock(0, 0, 0, 0)],
    FloquetExpansions.CKOutputConstraint{Int}[],
    [FloquetExpansions.CKOutputConstraint(1, 0, harmonic)],
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
ck_normalization_channel_pair(left, right) = kron(conj.(right), left)

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

  @test FloquetExpansions.ck_output_pairing_period_terms(normalized.coefficients[1]) ==
    Dict(0 => identity_component)
  @test all(isempty(coefficient.terms) for coefficient in normalized.coefficients[2:end])

  expected_s1 = ck_normalization_polynomial(
    0 => (-1 // 2) * m1_zero, 1 => (-1 // 2) * m1_one
  )
  @test normalization.coefficients[2].terms == expected_s1.terms

  for coefficient in normalization.coefficients, matrix in values(coefficient.terms)
    @test matrix == adjoint(matrix)
  end
end

@testset "phase-resolved CK normalization convolves Floquet grades" begin
  identity_component = ck_normalization_identity()
  zero_component = zero(identity_component)
  m_plus = CKNormalizationExact[0 1 + ck_normalization_im; 0 0]
  m_minus = adjoint(m_plus)

  metric = ck_normalization_series(
    ck_normalization_polynomial(0 => identity_component),
    ck_normalization_graded_polynomial((1, 0) => m_plus, (-1, 0) => m_minus),
  )
  normalization = FloquetExpansions.ck_output_metric_inverse_sqrt_series(
    metric, 2, identity_component
  )
  normalized = FloquetExpansions.ck_output_normalized_metric_series(
    metric, normalization, 2, zero_component
  )

  @test FloquetExpansions.ck_output_pairing_coefficient(
    normalization.coefficients[2], 1, 0
  ) == (-1 // 2) * m_plus
  @test FloquetExpansions.ck_output_pairing_coefficient(
    normalization.coefficients[2], -1, 0
  ) == (-1 // 2) * m_minus
  @test adjoint(
    FloquetExpansions.ck_output_pairing_coefficient(normalization.coefficients[2], 1, 0)
  ) == FloquetExpansions.ck_output_pairing_coefficient(
    normalization.coefficients[2], -1, 0
  )
  @test all(isempty(coefficient.terms) for coefficient in normalized.coefficients[2:end])

  vacuum = ck_normalization_vacuum_key()
  amplitude = ck_normalization_amplitude(0, vacuum => identity_component)
  normalized_amplitudes = FloquetExpansions.ck_output_right_normalize_series(
    [amplitude], normalization, 1, zero_component
  )
  phase_support = sort!(
    Int[key.phase_harmonic for key in keys(normalized_amplitudes[2].coefficients[1].terms)]
  )
  @test phase_support == [-1, 1]
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

  @test FloquetExpansions.ck_output_pairing_period_terms(normalization.coefficients[2]) ==
    Dict(3 => (-1 // 2) * m1)
  @test !iszero(
    FloquetExpansions.ck_output_pairing_coefficient(normalization.coefficients[3], 0, 6)
  )
  @test !iszero(
    FloquetExpansions.ck_output_pairing_coefficient(normalization.coefficients[4], 0, 9)
  )
  @test FloquetExpansions.ck_output_pairing_period_terms(normalized.coefficients[1]) ==
    Dict(0 => identity_component)
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
  @test FloquetExpansions.ck_output_pairing_period_terms(
    paired_normalized.coefficients[1]
  ) == Dict(0 => identity_component)
  @test all(
    isempty(coefficient.terms) for coefficient in paired_normalized.coefficients[2:end]
  )
  @test all(
    key.output_channels == [1] for coefficient in normalized_amplitudes for
    kernel in coefficient.coefficients for
    key in keys(kernel.terms) if !isempty(key.output_channels)
  )
end

@testset "singular physical Gram normalizes without Gram inversion" begin
  identity_component = ck_normalization_identity()
  zero_component = zero(identity_component)
  zero_superoperator = zeros(CKNormalizationExact, 4, 4)
  identity_superoperator = Matrix{CKNormalizationExact}(I, 4, 4)
  vacuum = ck_normalization_vacuum_key()
  phase = ck_normalization_scalar_resolvent_key(1)

  keys = [vacuum, phase]
  gram_right_left = CKNormalizationExact[
    FloquetExpansions.ck_output_overlap_coefficient(
      FloquetExpansions.ck_output_time_overlap(left, right, ck_normalization_im), 0
    ) for left in keys, right in keys
  ]
  @test gram_right_left ==
    CKNormalizationExact[1 -ck_normalization_im; ck_normalization_im 1]
  @test gram_right_left[1, 1] * gram_right_left[2, 2] ==
    gram_right_left[1, 2] * gram_right_left[2, 1]

  amplitude = ck_normalization_amplitude(
    0, vacuum => identity_component, phase => -ck_normalization_im * identity_component
  )
  raw_metric = FloquetExpansions.ck_output_metric_pairing(
    amplitude, amplitude, ck_normalization_im, zero_component, ck_normalization_metric_pair
  )
  raw_channel = FloquetExpansions.ck_output_channel_pairing(
    amplitude,
    amplitude,
    ck_normalization_im,
    zero_superoperator,
    ck_normalization_channel_pair,
  )
  @test FloquetExpansions.ck_output_pairing_period_terms(raw_metric) ==
    Dict(0 => 4 * identity_component)
  @test FloquetExpansions.ck_output_pairing_period_terms(raw_channel) ==
    Dict(0 => 4 * identity_superoperator)

  half_identity = (1 // 2) * identity_component
  normalization = FloquetExpansions.CKOutputPairingSeries([
    FloquetExpansions.CKOutputPairingPolynomial(Dict(0 => half_identity), zero_component)
  ])
  normalized_amplitude = only(
    FloquetExpansions.ck_output_right_normalize_series(
      [amplitude], normalization, 0, zero_component
    ),
  )
  normalized_metric = FloquetExpansions.ck_output_metric_pairing(
    normalized_amplitude,
    normalized_amplitude,
    ck_normalization_im,
    zero_component,
    ck_normalization_metric_pair,
  )
  normalized_channel = FloquetExpansions.ck_output_channel_pairing(
    normalized_amplitude,
    normalized_amplitude,
    ck_normalization_im,
    zero_superoperator,
    ck_normalization_channel_pair,
  )
  @test FloquetExpansions.ck_output_pairing_period_terms(normalized_metric) ==
    Dict(0 => identity_component)
  @test FloquetExpansions.ck_output_pairing_period_terms(normalized_channel) ==
    Dict(0 => identity_superoperator)
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

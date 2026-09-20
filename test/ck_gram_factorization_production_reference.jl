using Test
using FloquetExpansions
using LinearAlgebra: I, kron

const CKFactorReferenceExact = Complex{Rational{Int}}
const ck_factor_reference_im = CKFactorReferenceExact(0 // 1, 1 // 1)

function ck_factor_reference_vacuum_key()
  return FloquetExpansions.CKOutputKernelKey(
    Int[],
    FloquetExpansions.CKOutputBlock[],
    FloquetExpansions.CKOutputConstraint{Int}[],
    FloquetExpansions.CKOutputConstraint{Int}[],
  )
end

function ck_factor_reference_phase_key()
  return FloquetExpansions.CKOutputKernelKey(
    Int[],
    [FloquetExpansions.CKOutputBlock(0, 0, 0, 0)],
    FloquetExpansions.CKOutputConstraint{Int}[],
    [FloquetExpansions.CKOutputConstraint(1, 0, 1)],
  )
end

function ck_factor_reference_amplitude(terms...)
  zero_component = zeros(CKFactorReferenceExact, 2, 2)
  kernel = FloquetExpansions.CKOutputKernel(
    Dict{FloquetExpansions.CKOutputKernelKey{Int},Matrix{CKFactorReferenceExact}}(terms),
    zero_component,
  )
  return FloquetExpansions.CKOutputPeriodPolynomial([kernel])
end

ck_factor_reference_metric_pair(left, right) = adjoint(left) * right
ck_factor_reference_channel_pair(left, right) = kron(conj.(right), left)

@testset "production reverse-bra pairing equals independent Kraus factorization" begin
  vacuum = ck_factor_reference_vacuum_key()
  phase = ck_factor_reference_phase_key()
  keys = [vacuum, phase]

  gram_right_left = CKFactorReferenceExact[
    FloquetExpansions.ck_output_overlap_coefficient(
      FloquetExpansions.ck_output_time_overlap(left, right, ck_factor_reference_im), 0
    ) for left in keys, right in keys
  ]
  @test gram_right_left ==
    CKFactorReferenceExact[1 -ck_factor_reference_im; ck_factor_reference_im 1]

  factor = reshape(CKFactorReferenceExact[1, ck_factor_reference_im], 1, 2)
  gram = adjoint(factor) * factor
  @test transpose(gram) == gram_right_left

  first = CKFactorReferenceExact[1 1; 0 -1]
  second = CKFactorReferenceExact[0 1 + ck_factor_reference_im; 1 -1]
  amplitude = ck_factor_reference_amplitude(vacuum => first, phase => second)

  zero_component = zeros(CKFactorReferenceExact, 2, 2)
  zero_superoperator = zeros(CKFactorReferenceExact, 4, 4)
  production_metric = FloquetExpansions.ck_output_metric_pairing(
    amplitude,
    amplitude,
    ck_factor_reference_im,
    zero_component,
    ck_factor_reference_metric_pair,
  )
  production_channel = FloquetExpansions.ck_output_channel_pairing(
    amplitude,
    amplitude,
    ck_factor_reference_im,
    zero_superoperator,
    ck_factor_reference_channel_pair,
  )

  kraus = first + ck_factor_reference_im * second
  @test FloquetExpansions.ck_output_pairing_period_terms(production_metric) ==
    Dict(0 => adjoint(kraus) * kraus)
  @test FloquetExpansions.ck_output_pairing_period_terms(production_channel) ==
    Dict(0 => kron(conj.(kraus), kraus))
  @test FloquetExpansions.ck_output_pairing_phase_support(production_metric) == [0]
  @test FloquetExpansions.ck_output_pairing_phase_support(production_channel) == [0]
end

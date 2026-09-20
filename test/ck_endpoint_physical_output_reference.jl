using Test
using FloquetExpansions
using LinearAlgebra: I

const CKEndpointOutputExact = Complex{Rational{Int}}
const ck_endpoint_output_im = CKEndpointOutputExact(0 // 1, 1 // 1)

function ck_endpoint_output_zero_matrix()
  return zeros(CKEndpointOutputExact, 2, 2)
end

function ck_endpoint_output_identity_matrix()
  return Matrix{CKEndpointOutputExact}(I, 2, 2)
end

function ck_endpoint_output_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_endpoint_output_im / mismatch
end

function ck_endpoint_output_period_coefficient(values, degree::Int, zero_component)
  degree >= 0 || throw(ArgumentError("period degree must be nonnegative"))
  degree + 1 <= length(values) || return zero_component
  return values[degree + 1]
end

function ck_endpoint_output_time_coefficients(polynomial, output_channels, sidebands)
  return [
    begin
      result = coefficient.zero_component
      for (key, value) in coefficient.terms
        key.output_channels == output_channels || continue
        realization = FloquetExpansions.ck_output_time_realization(
          key, ck_endpoint_output_im
        )
        weight = FloquetExpansions.ck_output_time_fourier_weight(
          realization, sidebands, ck_endpoint_output_im
        )
        result += weight * value
      end
      result
    end for coefficient in polynomial.coefficients
  ]
end

function ck_endpoint_output_select_channels(polynomial, output_channels)
  coefficients = [
    FloquetExpansions.CKOutputKernel(
      Dict(
        key => value for
        (key, value) in coefficient.terms if key.output_channels == output_channels
      ),
      coefficient.zero_component,
    ) for coefficient in polynomial.coefficients
  ]
  return FloquetExpansions.CKOutputPeriodPolynomial(coefficients)
end

ck_endpoint_output_metric_pair(left, right) = adjoint(left) * right

@testset "endpoint reconstruction integrates one drift around one physical jump" begin
  zero_component = ck_endpoint_output_zero_matrix()
  identity_component = ck_endpoint_output_identity_matrix()

  jump_harmonic = 0
  drift_harmonic = 2
  jump_value = CKEndpointOutputExact[0 1; 0 0]
  drift_value = CKEndpointOutputExact[0 0; 1 0]

  A1 = FloquetExpansions.ck_kernel_generator(
    Dict(FloquetExpansions.ck_jump_vertex(1, jump_harmonic) => jump_value), zero_component
  )
  A2 = FloquetExpansions.ck_kernel_generator(
    Dict(FloquetExpansions.ck_drift_vertex(drift_harmonic) => drift_value), zero_component
  )
  identity_state = FloquetExpansions.ck_kernel_identity(
    0, identity_component, zero_component
  )
  zero_state = FloquetExpansions.ck_kernel_zero(0, zero_component)
  operations = FloquetExpansions.BlochOrderOperations(
    FloquetExpansions.ck_kernel_product,
    FloquetExpansions.ck_kernel_project_model,
    FloquetExpansions.ck_kernel_solve_complement,
  )

  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1, A2], 3, identity_state, zero_state, operations
  )
  reconstruction = FloquetExpansions.evaluate_ck_period_amplitude(
    recurrence.effective, recurrence.wave, 3, identity_state, zero_state
  )
  one_output = reconstruction.amplitude[4]
  finite_output = FloquetExpansions.ck_output_kernel(one_output)

  @test all(
    FloquetExpansions.ck_output_kernel_complete(key) for
    coefficient in finite_output.coefficients for key in keys(coefficient.terms)
  )
  @test all(
    key.phase_harmonic == drift_harmonic + jump_harmonic for
    coefficient in finite_output.coefficients for
    key in keys(coefficient.terms) if key.output_channels == [1]
  )

  before = jump_value * drift_value
  after = drift_value * jump_value
  expected_at_jump =
    (-ck_endpoint_output_im / drift_harmonic) * before +
    (ck_endpoint_output_im / drift_harmonic) * after
  expected_shifted = -expected_at_jump

  at_jump = FloquetExpansions.ck_period_ordered_sideband_coefficients(
    one_output, [1], [jump_harmonic]; inverse_weight=ck_endpoint_output_inverse_weight
  )
  at_shifted = FloquetExpansions.ck_period_ordered_sideband_coefficients(
    one_output,
    [1],
    [jump_harmonic + drift_harmonic];
    inverse_weight=ck_endpoint_output_inverse_weight,
  )
  outside = FloquetExpansions.ck_period_ordered_sideband_coefficients(
    one_output, [1], [jump_harmonic - 1]; inverse_weight=ck_endpoint_output_inverse_weight
  )
  finite_at_jump = FloquetExpansions.ck_output_ordered_sideband_coefficients(
    finite_output, [1], [jump_harmonic]; inverse_weight=ck_endpoint_output_inverse_weight
  )
  finite_at_shifted = FloquetExpansions.ck_output_ordered_sideband_coefficients(
    finite_output,
    [1],
    [jump_harmonic + drift_harmonic];
    inverse_weight=ck_endpoint_output_inverse_weight,
  )
  finite_outside = FloquetExpansions.ck_output_ordered_sideband_coefficients(
    finite_output,
    [1],
    [jump_harmonic - 1];
    inverse_weight=ck_endpoint_output_inverse_weight,
  )
  time_at_jump = ck_endpoint_output_time_coefficients(finite_output, [1], [jump_harmonic])
  time_at_shifted = ck_endpoint_output_time_coefficients(
    finite_output, [1], [jump_harmonic + drift_harmonic]
  )
  time_outside = ck_endpoint_output_time_coefficients(
    finite_output, [1], [jump_harmonic - 1]
  )

  @test finite_at_jump == at_jump
  @test finite_at_shifted == at_shifted
  @test finite_outside == outside
  @test time_at_jump == finite_at_jump
  @test time_at_shifted == finite_at_shifted
  @test time_outside == finite_outside
  @test ck_endpoint_output_period_coefficient(at_jump, 1, zero_component) ==
    expected_at_jump
  @test ck_endpoint_output_period_coefficient(at_shifted, 1, zero_component) ==
    expected_shifted
  @test all(iszero, ck_endpoint_output_period_coefficient(at_jump, 0, zero_component))
  @test all(iszero, ck_endpoint_output_period_coefficient(at_shifted, 0, zero_component))
  @test all(all(iszero, coefficient) for coefficient in outside)

  physical_one_output = ck_endpoint_output_select_channels(finite_output, [1])
  metric = FloquetExpansions.ck_output_metric_pairing(
    physical_one_output,
    physical_one_output,
    ck_endpoint_output_im,
    zero_component,
    ck_endpoint_output_metric_pair,
  )
  expected_metric = (2 // drift_harmonic^2) * identity_component
  @test FloquetExpansions.ck_output_pairing_period_terms(metric) == Dict(1 => expected_metric)
  @test FloquetExpansions.ck_output_pairing_phase_support(metric) == [0]
  @test !FloquetExpansions.ck_output_pairing_has_negative_power(metric)

  # The two matrix orderings isolate the raw-Dyson drift-before and drift-after
  # wavepackets from #197 without reintroducing their drift integration times.
  @test before == CKEndpointOutputExact[1 0; 0 0]
  @test after == CKEndpointOutputExact[0 0; 0 1]
  expected_matrix = CKEndpointOutputExact[
    -ck_endpoint_output_im / drift_harmonic 0
    0 ck_endpoint_output_im / drift_harmonic
  ]
  @test expected_at_jump == expected_matrix
end

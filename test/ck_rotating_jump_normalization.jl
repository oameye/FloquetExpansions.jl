using Test
using FloquetExpansions
using LinearAlgebra: I, kron

const CKRotatingNormalizationExact = Complex{Rational{Int}}
const ck_rotating_normalization_im = CKRotatingNormalizationExact(0 // 1, 1 // 1)

function ck_rotating_normalization_fixture()
  sigma_x = CKRotatingNormalizationExact[0 1; 1 0]
  sigma_y = CKRotatingNormalizationExact[
    0 -ck_rotating_normalization_im
    ck_rotating_normalization_im 0
  ]
  sigma_z = CKRotatingNormalizationExact[1 0; 0 -1]
  sigma_minus = (1 // 2) * (sigma_x - ck_rotating_normalization_im * sigma_y)
  sigma_plus = (1 // 2) * (sigma_x + ck_rotating_normalization_im * sigma_y)
  identity_component = Matrix{CKRotatingNormalizationExact}(I, 2, 2)
  zero_component = zero(identity_component)

  A1 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_jump_vertex(1, 0) => sigma_z,
      FloquetExpansions.ck_jump_vertex(1, 1) => ck_rotating_normalization_im * sigma_minus,
      FloquetExpansions.ck_jump_vertex(1, -1) => -ck_rotating_normalization_im * sigma_plus,
    ),
    zero_component,
  )
  # For this rotating jump L(t)^dagger L(t) = 2I, hence Q = -I after
  # factoring the fixed rate gamma.
  A2 = FloquetExpansions.ck_kernel_generator(
    Dict(FloquetExpansions.ck_drift_vertex(0) => -identity_component), zero_component
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
  return (;
    A1,
    A2,
    identity_state,
    zero_state,
    operations,
    identity_component,
    zero_component,
    sigma_z,
  )
end

function ck_rotating_normalization_amplitudes(order::Int)
  fixture = ck_rotating_normalization_fixture()
  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    order,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  reconstruction = FloquetExpansions.evaluate_ck_period_amplitude(
    recurrence.effective, recurrence.wave, order, fixture.identity_state, fixture.zero_state
  )
  amplitudes = [
    FloquetExpansions.ck_output_kernel(amplitude) for amplitude in reconstruction.amplitude
  ]
  return fixture, amplitudes
end

ck_rotating_normalization_metric_pair(left, right) = adjoint(left) * right
ck_rotating_normalization_channel_pair(left, right) = kron(conj.(right), left)

function ck_rotating_normalization_hamiltonian_super(H)
  identity_component = Matrix{CKRotatingNormalizationExact}(I, 2, 2)
  return -ck_rotating_normalization_im *
         (kron(identity_component, H) - kron(transpose(H), identity_component))
end

@testset "rotating-jump physical CK amplitudes are TP through the retained order" begin
  fixture, amplitudes = ck_rotating_normalization_amplitudes(4)
  metric = FloquetExpansions.ck_output_metric_series(
    amplitudes,
    4,
    ck_rotating_normalization_im,
    fixture.zero_component,
    ck_rotating_normalization_metric_pair,
  )

  @test FloquetExpansions.ck_output_pairing_period_terms(metric.coefficients[1]) ==
    Dict(0 => fixture.identity_component)
  @test FloquetExpansions.ck_output_pairing_phase_support(metric.coefficients[1]) == [0]
  @test all(isempty(coefficient.terms) for coefficient in metric.coefficients[2:end])
  @test all(
    !FloquetExpansions.ck_output_pairing_has_negative_power(coefficient) for
    coefficient in metric.coefficients
  )
end

@testset "new finite-output pairing reproduces #204 rotating-jump phase channel" begin
  fixture, amplitudes = ck_rotating_normalization_amplitudes(4)
  zero_superoperator = zeros(CKRotatingNormalizationExact, 4, 4)
  channel = FloquetExpansions.ck_output_channel_series(
    amplitudes,
    4,
    ck_rotating_normalization_im,
    zero_superoperator,
    ck_rotating_normalization_channel_pair,
  )

  raw_second_phase_support = sort!(unique([
    key.phase_harmonic for coefficient in amplitudes[5].coefficients for
    key in keys(coefficient.terms)
  ]))
  first_channel = channel.coefficients[3]
  second_channel = channel.coefficients[5]
  @test raw_second_phase_support == [-2, -1, 0, 1, 2]
  @test FloquetExpansions.ck_output_pairing_phase_support(second_channel) == [-1, 0, 1]
  @test iszero(FloquetExpansions.ck_output_pairing_coefficient(second_channel, -2, 1))
  @test iszero(FloquetExpansions.ck_output_pairing_coefficient(second_channel, 2, 1))
  @test !iszero(FloquetExpansions.ck_output_pairing_coefficient(second_channel, -1, 1))
  @test !iszero(FloquetExpansions.ck_output_pairing_coefficient(second_channel, 1, 1))

  first_square = FloquetExpansions.ck_output_pairing_product(
    first_channel, first_channel, zero_superoperator
  )
  log_second = FloquetExpansions.CKOutputPairingPolynomial(
    copy(second_channel.terms), zero_superoperator
  )
  FloquetExpansions.ck_output_pairing_add!(
    log_second,
    FloquetExpansions.ck_output_pairing_scale(
      -1 // 2, first_square, zero_superoperator
    ),
  )

  @test iszero(FloquetExpansions.ck_output_pairing_coefficient(log_second, 0, 2))
  expected = ck_rotating_normalization_hamiltonian_super(-(5 // 4) * fixture.sigma_z)
  @test FloquetExpansions.ck_output_pairing_coefficient(log_second, 0, 1) == expected
end

@testset "right normalization leaves retained rotating-jump physics unchanged" begin
  fixture, amplitudes = ck_rotating_normalization_amplitudes(4)
  metric = FloquetExpansions.ck_output_metric_series(
    amplitudes,
    4,
    ck_rotating_normalization_im,
    fixture.zero_component,
    ck_rotating_normalization_metric_pair,
  )
  normalization = FloquetExpansions.ck_output_metric_inverse_sqrt_series(
    metric, 4, fixture.identity_component
  )

  @test FloquetExpansions.ck_output_pairing_period_terms(normalization.coefficients[1]) ==
    Dict(0 => fixture.identity_component)
  @test FloquetExpansions.ck_output_pairing_phase_support(normalization.coefficients[1]) == [0]
  @test all(isempty(coefficient.terms) for coefficient in normalization.coefficients[2:end])

  normalized_amplitudes = FloquetExpansions.ck_output_right_normalize_series(
    amplitudes, normalization, 4, fixture.zero_component
  )
  normalized_metric = FloquetExpansions.ck_output_metric_series(
    normalized_amplitudes,
    4,
    ck_rotating_normalization_im,
    fixture.zero_component,
    ck_rotating_normalization_metric_pair,
  )
  @test normalized_metric.coefficients == metric.coefficients

  zero_superoperator = zeros(CKRotatingNormalizationExact, 4, 4)
  raw_channel = FloquetExpansions.ck_output_channel_series(
    amplitudes,
    4,
    ck_rotating_normalization_im,
    zero_superoperator,
    ck_rotating_normalization_channel_pair,
  )
  normalized_channel = FloquetExpansions.ck_output_channel_series(
    normalized_amplitudes,
    4,
    ck_rotating_normalization_im,
    zero_superoperator,
    ck_rotating_normalization_channel_pair,
  )
  @test normalized_channel.coefficients == raw_channel.coefficients
end

@testset "CK physical reconstruction is prefix-consistent through normalization" begin
  fixture3, amplitudes3 = ck_rotating_normalization_amplitudes(3)
  fixture4, amplitudes4 = ck_rotating_normalization_amplitudes(4)
  zero_superoperator = zeros(CKRotatingNormalizationExact, 4, 4)

  raw_channel3 = FloquetExpansions.ck_output_channel_series(
    amplitudes3,
    3,
    ck_rotating_normalization_im,
    zero_superoperator,
    ck_rotating_normalization_channel_pair,
  )
  raw_channel4_prefix = FloquetExpansions.ck_output_channel_series(
    amplitudes4,
    3,
    ck_rotating_normalization_im,
    zero_superoperator,
    ck_rotating_normalization_channel_pair,
  )
  @test raw_channel3.coefficients == raw_channel4_prefix.coefficients

  metric3 = FloquetExpansions.ck_output_metric_series(
    amplitudes3,
    3,
    ck_rotating_normalization_im,
    fixture3.zero_component,
    ck_rotating_normalization_metric_pair,
  )
  metric4_prefix = FloquetExpansions.ck_output_metric_series(
    amplitudes4,
    3,
    ck_rotating_normalization_im,
    fixture4.zero_component,
    ck_rotating_normalization_metric_pair,
  )
  @test metric3.coefficients == metric4_prefix.coefficients

  normalization3 = FloquetExpansions.ck_output_metric_inverse_sqrt_series(
    metric3, 3, fixture3.identity_component
  )
  normalization4_prefix = FloquetExpansions.ck_output_metric_inverse_sqrt_series(
    metric4_prefix, 3, fixture4.identity_component
  )
  @test normalization3.coefficients == normalization4_prefix.coefficients

  normalized3 = FloquetExpansions.ck_output_right_normalize_series(
    amplitudes3, normalization3, 3, fixture3.zero_component
  )
  normalized4_prefix = FloquetExpansions.ck_output_right_normalize_series(
    amplitudes4, normalization4_prefix, 3, fixture4.zero_component
  )
  normalized_channel3 = FloquetExpansions.ck_output_channel_series(
    normalized3,
    3,
    ck_rotating_normalization_im,
    zero_superoperator,
    ck_rotating_normalization_channel_pair,
  )
  normalized_channel4_prefix = FloquetExpansions.ck_output_channel_series(
    normalized4_prefix,
    3,
    ck_rotating_normalization_im,
    zero_superoperator,
    ck_rotating_normalization_channel_pair,
  )
  @test normalized_channel3.coefficients == normalized_channel4_prefix.coefficients
end

using Test
using FloquetExpansions
using LinearAlgebra: I, kron

const CKChannelGaugeExact = Complex{Rational{Int}}
const ck_channel_gauge_im = CKChannelGaugeExact(0 // 1, 1 // 1)

function ck_channel_gauge_fixture()
  H0 = CKChannelGaugeExact[2 1 + ck_channel_gauge_im; 1 - ck_channel_gauge_im -1]
  H1 = CKChannelGaugeExact[1 2 - ck_channel_gauge_im; -1 1 + ck_channel_gauge_im]
  H2 = CKChannelGaugeExact[ck_channel_gauge_im 1; 2 -1]
  hamiltonian = Dict(
    0 => H0,
    1 => H1,
    -1 => Matrix(adjoint(H1)),
    2 => H2,
    -2 => Matrix(adjoint(H2)),
  )
  jumps = Dict(
    -1 => CKChannelGaugeExact[1 0; 2 ck_channel_gauge_im],
    0 => CKChannelGaugeExact[0 1; -1 2],
    2 => CKChannelGaugeExact[1 - ck_channel_gauge_im 2; 0 -1],
  )

  zero_component = zeros(CKChannelGaugeExact, 2, 2)
  identity_component = Matrix{CKChannelGaugeExact}(I, 2, 2)
  A1 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_jump_vertex(1, harmonic) => value for (harmonic, value) in jumps
    ),
    zero_component,
  )
  A2 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_drift_vertex(harmonic) => -ck_channel_gauge_im * value for
      (harmonic, value) in hamiltonian
    ),
    zero_component,
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
  return (; hamiltonian, jumps, A1, A2, identity_state, zero_state, operations)
end

function ck_channel_gauge_hd_period_amplitudes(hd, identity_state, zero_state, order::Int)
  zero_endpoint = FloquetExpansions.ck_endpoint_kernel(zero_state)
  identity_endpoint = FloquetExpansions.ck_endpoint_kernel(identity_state)
  zero_period = FloquetExpansions.ck_period_zero(zero_endpoint)
  identity_period = FloquetExpansions.ck_period_constant(identity_endpoint)

  generator_series = [zero_period for _ in 0:order]
  for n in 1:order
    generator_series[n + 1] = FloquetExpansions.ck_period_constant(hd.generator[n])
  end
  endpoint_wave = FloquetExpansions.ck_zero_constant_series_exponential(
    generator_series, order, identity_period, zero_period
  )
  endpoint_wave_inverse = FloquetExpansions.ck_unit_series_inverse(
    endpoint_wave, order, identity_period, zero_period
  )

  slow_generator = [zero_period for _ in 0:order]
  for n in 1:order
    slow_generator[n + 1] = FloquetExpansions.ck_period_monomial(hd.effective[n], 1)
  end
  slow_propagator = FloquetExpansions.ck_zero_constant_series_exponential(
    slow_generator, order, identity_period, zero_period
  )
  dressed_left = FloquetExpansions.ck_truncated_series_product(
    endpoint_wave, slow_propagator, order, zero_period
  )
  return FloquetExpansions.ck_truncated_series_product(
    dressed_left, endpoint_wave_inverse, order, zero_period
  )
end

function ck_channel_gauge_output_amplitudes(amplitudes)
  return [FloquetExpansions.ck_output_kernel(amplitude) for amplitude in amplitudes]
end

ck_channel_gauge_pair(left, right) = kron(conj(right), left)

function ck_channel_gauge_channel_series(amplitudes, order::Int)
  output_amplitudes = ck_channel_gauge_output_amplitudes(amplitudes)
  zero_channel = zeros(CKChannelGaugeExact, 4, 4)
  return FloquetExpansions.ck_output_channel_series(
    output_amplitudes,
    order,
    ck_channel_gauge_im,
    zero_channel,
    ck_channel_gauge_pair,
  )
end

@testset "canonical BF and direct HD agree after exact physical channel pasting" begin
  fixture = ck_channel_gauge_fixture()
  order = 5

  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    order,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  canonical = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence.effective,
    recurrence.wave,
    order,
    fixture.identity_state,
    fixture.zero_state,
  )
  bf_period = FloquetExpansions.evaluate_ck_period_amplitude(
    canonical.effective,
    canonical.wave,
    order,
    fixture.identity_state,
    fixture.zero_state,
  )

  hd = FloquetExpansions.evaluate_ck_hori_deprit(
    [fixture.A1, fixture.A2], order, fixture.zero_state
  )
  hd_amplitudes = ck_channel_gauge_hd_period_amplitudes(
    hd, fixture.identity_state, fixture.zero_state, order
  )

  bf_channel = ck_channel_gauge_channel_series(bf_period.amplitude, order)
  hd_channel = ck_channel_gauge_channel_series(hd_amplitudes, order)

  @test bf_channel.coefficients == hd_channel.coefficients
  @test any(!isempty(coefficient.terms) for coefficient in bf_channel.coefficients[2:end])
end

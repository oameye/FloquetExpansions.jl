using Test
using FloquetExpansions
using LinearAlgebra: I, kron

const CKChannelGaugeExact = Complex{Rational{Int}}
const ck_channel_gauge_im = CKChannelGaugeExact(0 // 1, 1 // 1)

function ck_channel_gauge_fixture()
  H0 = CKChannelGaugeExact[1 1; 1 -1]
  H1 = CKChannelGaugeExact[1 1 + ck_channel_gauge_im; 2 -1]
  hamiltonian = Dict(0 => H0, 1 => H1, -1 => Matrix(adjoint(H1)))
  jumps = Dict(
    0 => CKChannelGaugeExact[0 1; 0 0],
    1 => CKChannelGaugeExact[1 0; 1 ck_channel_gauge_im],
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

function ck_channel_gauge_slow_amplitudes(effective, identity_state, zero_state, order::Int)
  zero_endpoint = FloquetExpansions.ck_endpoint_kernel(zero_state)
  identity_endpoint = FloquetExpansions.ck_endpoint_kernel(identity_state)
  zero_period = FloquetExpansions.ck_period_zero(zero_endpoint)
  identity_period = FloquetExpansions.ck_period_constant(identity_endpoint)
  slow_generator = [zero_period for _ in 0:order]
  for n in 1:order
    slow_generator[n + 1] = FloquetExpansions.ck_period_monomial(effective[n], 1)
  end
  return FloquetExpansions.ck_zero_constant_series_exponential(
    slow_generator, order, identity_period, zero_period
  )
end

function ck_channel_gauge_one_output(amplitude)
  output = FloquetExpansions.ck_output_kernel(amplitude)
  coefficients = typeof(output.coefficients)(undef, length(output.coefficients))
  for index in eachindex(output.coefficients)
    kernel = output.coefficients[index]
    terms = Dict(
      key => value for (key, value) in kernel.terms if length(key.output_channels) == 1
    )
    coefficients[index] = FloquetExpansions.CKOutputKernel(terms, kernel.zero_component)
  end
  return FloquetExpansions.CKOutputPeriodPolynomial(coefficients)
end

function ck_channel_gauge_one_output_series(amplitudes)
  return [ck_channel_gauge_one_output(amplitude) for amplitude in amplitudes]
end

ck_channel_gauge_pair(left, right) = kron(conj(right), left)

function ck_channel_gauge_channel_series(amplitudes, order::Int)
  zero_channel = zeros(CKChannelGaugeExact, 4, 4)
  return FloquetExpansions.ck_output_channel_series(
    ck_channel_gauge_one_output_series(amplitudes),
    order,
    ck_channel_gauge_im,
    zero_channel,
    ck_channel_gauge_pair,
  )
end

@testset "canonical BF and direct HD agree after one-dissipator physical pasting" begin
  fixture = ck_channel_gauge_fixture()
  amplitude_order = 5
  channel_order = 6

  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    amplitude_order,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  canonical = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence.effective,
    recurrence.wave,
    amplitude_order,
    fixture.identity_state,
    fixture.zero_state,
  )
  hd = FloquetExpansions.evaluate_ck_hori_deprit(
    [fixture.A1, fixture.A2], amplitude_order, fixture.zero_state
  )

  bf_slow = ck_channel_gauge_slow_amplitudes(
    canonical.effective, fixture.identity_state, fixture.zero_state, amplitude_order
  )
  hd_slow = ck_channel_gauge_slow_amplitudes(
    hd.effective, fixture.identity_state, fixture.zero_state, amplitude_order
  )
  bf_channel = ck_channel_gauge_channel_series(bf_slow, channel_order)
  hd_channel = ck_channel_gauge_channel_series(hd_slow, channel_order)

  @test bf_channel.coefficients == hd_channel.coefficients
  @test isempty(bf_channel.coefficients[1].terms)
  @test isempty(bf_channel.coefficients[2].terms)
  @test !isempty(bf_channel.coefficients[3].terms)
  @test !isempty(bf_channel.coefficients[7].terms)
end

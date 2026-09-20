using Test
using FloquetExpansions
using LinearAlgebra: I, kron

const CKChannelGaugeExact = Complex{Rational{Int}}
const ck_channel_gauge_im = CKChannelGaugeExact(0 // 1, 1 // 1)

include("helpers/ck_channel_gauge_reference.jl")

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
  return (; hamiltonian, jumps, A1, A2, identity_state, zero_state, zero_component, operations)
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

function ck_channel_gauge_output_number(amplitude, output_number::Int)
  output = FloquetExpansions.ck_output_kernel(amplitude)
  coefficients = typeof(output.coefficients)(undef, length(output.coefficients))
  for index in eachindex(output.coefficients)
    kernel = output.coefficients[index]
    terms = Dict(
      key => value for (key, value) in kernel.terms if
      length(key.output_channels) == output_number
    )
    coefficients[index] = FloquetExpansions.CKOutputKernel(terms, kernel.zero_component)
  end
  return FloquetExpansions.CKOutputPeriodPolynomial(coefficients)
end

function ck_channel_gauge_output_series(amplitudes, output_number::Int)
  return [
    ck_channel_gauge_output_number(amplitude, output_number) for amplitude in amplitudes
  ]
end

ck_channel_gauge_pair(left, right) = kron(conj(right), left)

function ck_channel_gauge_channel_series(amplitudes, order::Int, output_number::Int)
  zero_channel = zeros(CKChannelGaugeExact, 4, 4)
  return FloquetExpansions.ck_output_channel_series(
    ck_channel_gauge_output_series(amplitudes, output_number),
    order,
    ck_channel_gauge_im,
    zero_channel,
    ck_channel_gauge_pair,
  )
end

function ck_channel_gauge_canonical_effective(fixture, A1, A2, order::Int)
  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1, A2],
    order,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  return FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence.effective,
    recurrence.wave,
    order,
    fixture.identity_state,
    fixture.zero_state,
  ).effective
end

function ck_channel_gauge_scaled_drift(fixture, scale::Int)
  loss = ck_channel_gauge_loss_harmonics(fixture.jumps, fixture.zero_component)
  harmonics = union(keys(fixture.hamiltonian), keys(loss))
  terms = Dict{FloquetExpansions.CKPhysicalVertex{Int},Matrix{CKChannelGaugeExact}}()
  for harmonic in harmonics
    value =
      -ck_channel_gauge_im *
      get(fixture.hamiltonian, harmonic, fixture.zero_component) +
      scale * get(loss, harmonic, fixture.zero_component)
    all(iszero, value) && continue
    terms[FloquetExpansions.ck_drift_vertex(harmonic)] = value
  end
  return FloquetExpansions.ck_kernel_generator(terms, fixture.zero_component)
end

function ck_channel_gauge_zero_output_order6(fixture, scale::Int)
  amplitude_order = 6
  A2 = ck_channel_gauge_scaled_drift(fixture, scale)
  bf_effective = ck_channel_gauge_canonical_effective(
    fixture, fixture.zero_state, A2, amplitude_order
  )
  hd = FloquetExpansions.evaluate_ck_hori_deprit(
    [fixture.zero_state, A2], amplitude_order, fixture.zero_state
  )
  bf_slow = ck_channel_gauge_slow_amplitudes(
    bf_effective, fixture.identity_state, fixture.zero_state, amplitude_order
  )
  hd_slow = ck_channel_gauge_slow_amplitudes(
    hd.effective, fixture.identity_state, fixture.zero_state, amplitude_order
  )
  bf_channel = ck_channel_gauge_channel_series(bf_slow, amplitude_order, 0)
  hd_channel = ck_channel_gauge_channel_series(hd_slow, amplitude_order, 0)
  bf_value = FloquetExpansions.ck_output_pairing_coefficient(
    bf_channel.coefficients[7], 0, 1
  )
  hd_value = FloquetExpansions.ck_output_pairing_coefficient(
    hd_channel.coefficients[7], 0, 1
  )
  return bf_value, hd_value
end

function ck_channel_gauge_linear_nojump_order6(fixture)
  values = Dict{Int,Matrix{CKChannelGaugeExact}}()
  for scale in (-2, -1, 1, 2)
    bf_value, hd_value = ck_channel_gauge_zero_output_order6(fixture, scale)
    @test bf_value == hd_value
    values[scale] = bf_value
  end
  first_difference = (values[1] - values[-1]) / 2
  second_difference = (values[2] - values[-2]) / 4
  return (4 * first_difference - second_difference) / 3
end

@testset "canonical BF and direct HD agree after one-dissipator physical pasting" begin
  fixture = ck_channel_gauge_fixture()
  amplitude_order = 5
  channel_order = 6
  bf_effective = ck_channel_gauge_canonical_effective(
    fixture, fixture.A1, fixture.A2, amplitude_order
  )
  hd = FloquetExpansions.evaluate_ck_hori_deprit(
    [fixture.A1, fixture.A2], amplitude_order, fixture.zero_state
  )

  bf_slow = ck_channel_gauge_slow_amplitudes(
    bf_effective, fixture.identity_state, fixture.zero_state, amplitude_order
  )
  hd_slow = ck_channel_gauge_slow_amplitudes(
    hd.effective, fixture.identity_state, fixture.zero_state, amplitude_order
  )
  bf_channel = ck_channel_gauge_channel_series(bf_slow, channel_order, 1)
  hd_channel = ck_channel_gauge_channel_series(hd_slow, channel_order, 1)

  @test bf_channel.coefficients == hd_channel.coefficients
  @test isempty(bf_channel.coefficients[1].terms)
  @test isempty(bf_channel.coefficients[2].terms)
  @test !isempty(bf_channel.coefficients[3].terms)
  @test !isempty(bf_channel.coefficients[7].terms)
end

@testset "canonical CK obeys the #104 one-dissipator static similarity after pasting" begin
  fixture = ck_channel_gauge_fixture()
  amplitude_order = 5
  channel_order = 6
  bf_effective = ck_channel_gauge_canonical_effective(
    fixture, fixture.A1, fixture.A2, amplitude_order
  )
  bf_slow = ck_channel_gauge_slow_amplitudes(
    bf_effective, fixture.identity_state, fixture.zero_state, amplitude_order
  )
  jump_channel = ck_channel_gauge_channel_series(bf_slow, channel_order, 1)
  jump_second_order = FloquetExpansions.ck_output_pairing_coefficient(
    jump_channel.coefficients[7], 0, 1
  )
  nojump_second_order = ck_channel_gauge_linear_nojump_order6(fixture)
  canonical_second_order = jump_second_order + nojump_second_order

  cp_second_order = ck_channel_gauge_cp_second_order(
    fixture.hamiltonian, fixture.jumps, fixture.zero_component
  )
  static_correction, H0_map = ck_channel_gauge_static_correction(
    fixture.hamiltonian, fixture.jumps, fixture.zero_component
  )
  gauge_shift = ck_channel_gauge_commutator(static_correction, H0_map)

  @test !all(iszero, static_correction)
  @test !all(iszero, gauge_shift)
  @test canonical_second_order == cp_second_order + gauge_shift
  @test canonical_second_order != cp_second_order
end

using Test
using FloquetExpansions

const CKOutputPairingExact = Complex{Rational{Int}}
const ck_output_pairing_im = CKOutputPairingExact(0 // 1, 1 // 1)

function ck_output_pairing_model_key(
  cumulative_harmonics::Vector{Int}, phase_harmonic::Int=0
)
  output_number = length(cumulative_harmonics)
  constraints = FloquetExpansions.CKOutputConstraint{Int}[
    FloquetExpansions.CKOutputConstraint(1, output_stop, harmonic) for
    (output_stop, harmonic) in enumerate(cumulative_harmonics)
  ]
  return FloquetExpansions.CKOutputKernelKey(
    fill(1, output_number),
    [FloquetExpansions.CKOutputBlock(0, output_number, 0, output_number)],
    constraints,
    FloquetExpansions.CKOutputConstraint{Int}[],
    phase_harmonic,
  )
end

function ck_output_pairing_period_polynomial(
  period_power::Int,
  amplitudes::Vector{Pair{FloquetExpansions.CKOutputKernelKey{Int},CKOutputPairingExact}},
)
  zero_component = zero(CKOutputPairingExact)
  coefficients = FloquetExpansions.CKOutputKernel{Int,CKOutputPairingExact}[
    FloquetExpansions.CKOutputKernel(
      Dict{FloquetExpansions.CKOutputKernelKey{Int},CKOutputPairingExact}(), zero_component
    ) for _ in 0:period_power
  ]
  coefficients[period_power + 1] = FloquetExpansions.CKOutputKernel(
    Dict(amplitudes), zero_component
  )
  return FloquetExpansions.CKOutputPeriodPolynomial(coefficients)
end

ck_output_pairing_channel(left, right) = left * conj(right)
ck_output_pairing_metric(left, right) = conj(left) * right

@testset "inverse Fourier normalization restores the physical one-output norm" begin
  key = ck_output_pairing_model_key([3])
  amplitude = ck_output_pairing_period_polynomial(1, [key => one(CKOutputPairingExact)])
  channel = FloquetExpansions.ck_output_channel_pairing(
    amplitude,
    amplitude,
    ck_output_pairing_im,
    zero(CKOutputPairingExact),
    ck_output_pairing_channel,
  )

  @test FloquetExpansions.ck_output_pairing_period_terms(channel) ==
    Dict(1 => one(CKOutputPairingExact))
  @test !FloquetExpansions.ck_output_pairing_has_negative_power(channel)
end

@testset "one-drift dressed jump has the frozen physical norm and interference" begin
  drift_harmonic = 2
  jump_harmonic = 3
  jump = ck_output_pairing_model_key([jump_harmonic])
  shifted = ck_output_pairing_model_key([jump_harmonic + drift_harmonic])
  before = ck_output_pairing_period_polynomial(
    1,
    [
      jump => -ck_output_pairing_im / drift_harmonic,
      shifted => ck_output_pairing_im / drift_harmonic,
    ],
  )
  after = ck_output_pairing_period_polynomial(
    1,
    [
      jump => ck_output_pairing_im / drift_harmonic,
      shifted => -ck_output_pairing_im / drift_harmonic,
    ],
  )

  before_norm = FloquetExpansions.ck_output_channel_pairing(
    before,
    before,
    ck_output_pairing_im,
    zero(CKOutputPairingExact),
    ck_output_pairing_channel,
  )
  interference = FloquetExpansions.ck_output_channel_pairing(
    before,
    after,
    ck_output_pairing_im,
    zero(CKOutputPairingExact),
    ck_output_pairing_channel,
  )

  expected_norm = CKOutputPairingExact(2 // drift_harmonic^2)
  @test FloquetExpansions.ck_output_pairing_period_terms(before_norm) ==
    Dict(1 => expected_norm)
  @test FloquetExpansions.ck_output_pairing_period_terms(interference) ==
    Dict(1 => -expected_norm)
  @test !FloquetExpansions.ck_output_pairing_has_negative_power(before_norm)
  @test !FloquetExpansions.ck_output_pairing_has_negative_power(interference)
end

@testset "reverse-bra Gram orientation is distinct for channel and metric contractions" begin
  left_key = ck_output_pairing_model_key([1, 0])
  right_key = ck_output_pairing_model_key([0, 0])
  left = ck_output_pairing_period_polynomial(2, [left_key => one(CKOutputPairingExact)])
  right = ck_output_pairing_period_polynomial(2, [right_key => one(CKOutputPairingExact)])

  channel = FloquetExpansions.ck_output_channel_pairing(
    left, right, ck_output_pairing_im, zero(CKOutputPairingExact), ck_output_pairing_channel
  )
  metric = FloquetExpansions.ck_output_metric_pairing(
    left, right, ck_output_pairing_im, zero(CKOutputPairingExact), ck_output_pairing_metric
  )

  @test FloquetExpansions.ck_output_pairing_period_terms(channel) ==
    Dict(1 => ck_output_pairing_im)
  @test FloquetExpansions.ck_output_pairing_period_terms(metric) ==
    Dict(1 => -ck_output_pairing_im)
  @test !FloquetExpansions.ck_output_pairing_has_negative_power(channel)
  @test !FloquetExpansions.ck_output_pairing_has_negative_power(metric)
end

@testset "channel and metric carry opposite Floquet phase grades" begin
  left_key = ck_output_pairing_model_key([1, 0], 2)
  right_key = ck_output_pairing_model_key([0, 0], -1)
  left = ck_output_pairing_period_polynomial(2, [left_key => one(CKOutputPairingExact)])
  right = ck_output_pairing_period_polynomial(2, [right_key => one(CKOutputPairingExact)])

  channel = FloquetExpansions.ck_output_channel_pairing(
    left, right, ck_output_pairing_im, zero(CKOutputPairingExact), ck_output_pairing_channel
  )
  metric = FloquetExpansions.ck_output_metric_pairing(
    left, right, ck_output_pairing_im, zero(CKOutputPairingExact), ck_output_pairing_metric
  )

  @test FloquetExpansions.ck_output_pairing_phase_support(channel) == [3]
  @test FloquetExpansions.ck_output_pairing_phase_support(metric) == [-3]
  @test FloquetExpansions.ck_output_pairing_coefficient(channel, 3, 1) == ck_output_pairing_im
  @test FloquetExpansions.ck_output_pairing_coefficient(metric, -3, 1) == -ck_output_pairing_im
end

@testset "system factors follow the authoritative #204 conjugation convention" begin
  left_key = ck_output_pairing_model_key([1, 0])
  right_key = ck_output_pairing_model_key([0, 0])
  left_value = CKOutputPairingExact(2 // 1, 1 // 1)
  right_value = CKOutputPairingExact(3 // 1, 2 // 1)
  left = ck_output_pairing_period_polynomial(2, [left_key => left_value])
  right = ck_output_pairing_period_polynomial(2, [right_key => right_value])

  channel = FloquetExpansions.ck_output_channel_pairing(
    left, right, ck_output_pairing_im, zero(CKOutputPairingExact), ck_output_pairing_channel
  )
  metric = FloquetExpansions.ck_output_metric_pairing(
    left, right, ck_output_pairing_im, zero(CKOutputPairingExact), ck_output_pairing_metric
  )

  gram_rl = ck_output_pairing_im
  @test FloquetExpansions.ck_output_pairing_period_terms(channel) ==
    Dict(1 => gram_rl * left_value * conj(right_value))
  @test FloquetExpansions.ck_output_pairing_period_terms(metric) ==
    Dict(1 => conj(gram_rl) * conj(left_value) * right_value)
end

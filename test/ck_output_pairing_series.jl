using Test
using FloquetExpansions

const CKPairingSeriesExact = Complex{Rational{Int}}
const ck_pairing_series_im = CKPairingSeriesExact(0 // 1, 1 // 1)

function ck_pairing_series_vacuum_key()
  return FloquetExpansions.CKOutputKernelKey(
    Int[],
    FloquetExpansions.CKOutputBlock[],
    FloquetExpansions.CKOutputConstraint{Int}[],
    FloquetExpansions.CKOutputConstraint{Int}[],
  )
end

function ck_pairing_series_jump_key(harmonic::Int)
  return FloquetExpansions.CKOutputKernelKey(
    [1],
    [FloquetExpansions.CKOutputBlock(0, 1, 0, 1)],
    [FloquetExpansions.CKOutputConstraint(1, 1, harmonic)],
    FloquetExpansions.CKOutputConstraint{Int}[],
  )
end

function ck_pairing_series_period_polynomial(period_power::Int, key, value)
  zero_component = zero(CKPairingSeriesExact)
  coefficients = FloquetExpansions.CKOutputKernel{Int,CKPairingSeriesExact}[
    FloquetExpansions.CKOutputKernel(
      Dict{FloquetExpansions.CKOutputKernelKey{Int},CKPairingSeriesExact}(), zero_component
    ) for _ in 0:period_power
  ]
  coefficients[period_power + 1] = FloquetExpansions.CKOutputKernel(
    Dict(key => value), zero_component
  )
  return FloquetExpansions.CKOutputPeriodPolynomial(coefficients)
end

ck_pairing_series_channel(left, right) = left * conj(right)
ck_pairing_series_metric(left, right) = conj(left) * right

@testset "perturbative Stinespring pairing convolves amplitude orders exactly" begin
  vacuum = ck_pairing_series_period_polynomial(
    0, ck_pairing_series_vacuum_key(), one(CKPairingSeriesExact)
  )
  jump = ck_pairing_series_period_polynomial(
    1, ck_pairing_series_jump_key(0), one(CKPairingSeriesExact)
  )
  drift = ck_pairing_series_period_polynomial(
    1, ck_pairing_series_vacuum_key(), CKPairingSeriesExact(-1 // 2)
  )
  amplitudes = [vacuum, jump, drift]

  metric = FloquetExpansions.ck_output_metric_series(
    amplitudes,
    2,
    ck_pairing_series_im,
    zero(CKPairingSeriesExact),
    ck_pairing_series_metric,
  )
  channel = FloquetExpansions.ck_output_channel_series(
    amplitudes,
    2,
    ck_pairing_series_im,
    zero(CKPairingSeriesExact),
    ck_pairing_series_channel,
  )

  @test length(metric.coefficients) == 3
  @test metric.coefficients[1].terms == Dict(0 => one(CKPairingSeriesExact))
  @test isempty(metric.coefficients[2].terms)
  @test isempty(metric.coefficients[3].terms)
  @test channel.coefficients == metric.coefficients
end

@testset "static scalar GKSL cancellation is lost exactly beyond the retained triangle" begin
  vacuum = ck_pairing_series_period_polynomial(
    0, ck_pairing_series_vacuum_key(), one(CKPairingSeriesExact)
  )
  jump = ck_pairing_series_period_polynomial(
    1, ck_pairing_series_jump_key(0), one(CKPairingSeriesExact)
  )
  drift = ck_pairing_series_period_polynomial(
    1, ck_pairing_series_vacuum_key(), CKPairingSeriesExact(-1 // 2)
  )
  amplitudes = [vacuum, jump, drift]

  metric = FloquetExpansions.ck_output_metric_series(
    amplitudes,
    4,
    ck_pairing_series_im,
    zero(CKPairingSeriesExact),
    ck_pairing_series_metric,
  )

  @test isempty(metric.coefficients[2].terms)
  @test isempty(metric.coefficients[3].terms)
  @test isempty(metric.coefficients[4].terms)
  @test metric.coefficients[5].terms == Dict(2 => CKPairingSeriesExact(1 // 4))
end

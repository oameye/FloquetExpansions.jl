using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const FE = FloquetExpansions
const SQA = SecondQuantizedAlgebra

liouvillian_vanishes(L::Liouvillian) = iszero(FE.canonical_liouvillian(L))

@testset "native GKSL expansion represents a periodic scalar rate by its collapse amplitude" begin
  fock = FockSpace(:native_gksl_periodic_rate)
  a = Destroy(fock, :a)
  @variables ω_native_rate::Real t_native_rate::Real

  ω = ω_native_rate
  t = t_native_rate

  # The periodic rate Γ(t) = (1 + cos(ωt)/2)² of D[a] is the square of the collapse amplitude
  # C(t) = (1 + cos(ωt)/2) a, so Γ(t) D[a] = D[C(t)]. Squaring doubles the Fourier support: the
  # rate has the harmonics -2 to 2 and the amplitude only -1 to 1. Rate harmonics are not
  # amplitude harmonics and must not be transported as if they carried the same grading.
  amplitude_factor = 1 + (1 // 2) * cos(ω * t)
  rate_map = harmonics(amplitude_factor^2 * dissipator(a), ω, t)

  # The native GKSL expansion therefore consumes the amplitude itself. Without a coherent drive
  # it builds one harmonic jump channel per Fourier harmonic of that amplitude, and its generator
  # is the period average of the periodic-rate dissipator.
  native = FE.native_gksl_expansion(0 * a, ω, t, 1, (collapse(amplitude_factor * a),))
  @test sort!([channel.harmonic for channel in native.channels]) == [-1, 0, 1]
  @test liouvillian_vanishes(native.generator - time_average(rate_map))
end

@testset "native GKSL expansion: a new dissipative direction has an even rate onset" begin
  fock = FockSpace(:native_gksl_even_onset)
  a = Destroy(fock, :a)
  frame = DissipativeFrame(a^2, a)
  @variables ω_native_onset::Real t_native_onset::Real

  ω = ω_native_onset
  t = t_native_onset

  # A linear coherent drive dresses two-photon loss. The first-order micromotion generates an `a`
  # amplitude at grade one, so under integer amplitude grading the positive one-photon rate can
  # start only at grade two, once the transported amplitude is squared. At order two that grade
  # is the top of the finite square, which the native GKSL expansion keeps untruncated.
  #
  # The odd counterpart, a rate onset at grade one in a direction dark at grade zero, requires a
  # fractional amplitude onset. Gram positive completion reports it as `FractionalJumpOnset`,
  # certified in `cp_completion_validation.jl`.
  H = cos(ω * t) * (a + a')
  native = FE.native_gksl_expansion(H, ω, t, 2, (collapse(a^2),))
  rate(grade) =
    kossakowski(FE.native_dissipative_component(native.amplitudes, grade), frame)

  @test iszero(SQA.simplify(rate(0)[2, 2]))
  @test iszero(SQA.simplify(rate(1)[2, 2]))
  @test !iszero(SQA.simplify(rate(2)[2, 2]))
end

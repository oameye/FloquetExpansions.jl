using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const FE = FloquetExpansions
const SQA = SecondQuantizedAlgebra

cp_rate_liouvillian_zero(L::Liouvillian) = iszero(FE.canonical_liouvillian(L))

@testset "periodic scalar rates are represented by physical collapse amplitudes" begin
  fock = FockSpace(:cp_hfe_periodic_rate)
  a = Destroy(fock, :a)
  @variables ω_cp_rate::Real t_cp_rate::Real γ_cp_rate::Real

  ω = ω_cp_rate
  t = t_cp_rate
  γ = γ_cp_rate

  amplitude_factor = 1 + (1 // 2) * cos(ω * t)
  amplitude = amplitude_factor * a
  rate = amplitude_factor^2

  amplitude_harmonics = harmonics(amplitude, ω, t)
  collapse_map = harmonics(dissipator(amplitude), ω, t)
  rate_map = harmonics(rate * dissipator(a), ω, t)

  # The physical rate is the Gram square of the collapse amplitude. Squaring doubles the
  # Fourier support: rate harmonics are not amplitude harmonics and must not be transported as if
  # they carried the same perturbative grading.
  @test sort!(collect(keys(amplitude_harmonics))) == [-1, 0, 1]
  @test sort!(collect(keys(rate_map))) == [-2, -1, 0, 1, 2]
  @test sort!(collect(keys(collapse_map))) == [-2, -1, 0, 1, 2]
  for harmonic in -2:2
    @test cp_rate_liouvillian_zero(collapse_map[harmonic] - rate_map[harmonic])
  end

  # Native CP-HFE therefore consumes the square-root object itself. With no coherent drive, the
  # finite Kraus rows are exactly the Fourier components of that physical collapse amplitude.
  native = FE.cp_hfe_reconstruction(0 * a, ω, t, 1, (collapse(amplitude),))
  @test sort!([channel.harmonic for channel in native.channels]) == [-1, 0, 1]
  @test cp_rate_liouvillian_zero(native.generator - time_average(rate_map))

  # A separate static physical rate remains an external nonnegative channel weight and does not
  # alter the amplitude harmonic structure.
  weighted = FE.cp_hfe_reconstruction(0 * a, ω, t, 1, (jump(amplitude, γ),))
  @test only(weighted.amplitudes).seed.rate == γ
  @test cp_rate_liouvillian_zero(weighted.generator - γ * native.generator)

  # Passing the squared periodic rate as `jump(a, rate(t))` would lose its amplitude square root.
  # The native path must reject that representation rather than assigning rate coefficients to
  # amplitude grades.
  error = try
    FE.cp_hfe_reconstruction(0 * a, ω, t, 1, (jump(a, rate),))
    nothing
  catch caught
    caught
  end
  @test error isa ArgumentError
  @test occursin("explicit periodic collapse amplitude", sprint(showerror, error))
end

@testset "integer amplitude grading gives even onset for a new dissipative direction" begin
  fock = FockSpace(:cp_hfe_even_amplitude_onset)
  a = Destroy(fock, :a)
  frame = DissipativeFrame(a^2, a)
  @variables ω_cp_onset::Real t_cp_onset::Real

  ω = ω_cp_onset
  t = t_cp_onset

  # A linear coherent drive dresses two-photon loss. The first kick generates an `a` amplitude at
  # grade one; its positive one-photon rate can therefore start only at grade two after the Gram
  # square is formed.
  H = cos(ω * t) * (a + a')
  native = FE.cp_hfe_reconstruction(H, ω, t, 2, (collapse(a^2),))
  amplitudes = native.amplitudes

  rate0 = kossakowski(FE.cp_dissipative_component(amplitudes, 0), frame)
  rate1 = kossakowski(FE.cp_dissipative_component(amplitudes, 1), frame)
  rate2 = kossakowski(FE.cp_dissipative_component(amplitudes, 2), frame)

  @test !iszero(SQA.simplify(rate0[1, 1]))
  @test iszero(SQA.simplify(rate0[2, 2]))
  @test iszero(SQA.simplify(rate1[2, 2]))
  @test !iszero(SQA.simplify(rate2[2, 2]))
end

@testset "odd dark-sector rate onset requires a fractional amplitude onset" begin
  pauli = PauliSpace(:cp_hfe_fractional_onset_control)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  @variables ω_cp_fractional::Real

  ω = ω_cp_fractional
  coherent_rotation = hamiltonian_action(σx)
  cross_dissipator = dissipator(σy + σz) - dissipator(σy) - dissipator(σz)
  quadrature = im * (1 // 2) * cross_dissipator
  generator = PeriodicGenerator(
    Dict(
      0 => dissipator(σz),
      1 => coherent_rotation + quadrature,
      -1 => coherent_rotation - quadrature,
    ),
    ω,
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)

  error = try
    positive_completion(expansion, Gram(), frame)
    nothing
  catch caught
    caught
  end
  @test error isa FE.FractionalJumpOnset
  @test error.rate_order == 1
end

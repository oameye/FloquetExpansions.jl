using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const FE = FloquetExpansions
const SQA = SecondQuantizedAlgebra

native_validation_liouvillian_zero(L::Liouvillian) = iszero(FE.canonical_liouvillian(L))

function native_validation_retained_component(native, wd, grade::Int)
  coherent = effective_component(native.coherent, grade)
  dissipative = FE.native_dissipative_component(native.amplitudes, grade)
  return SQA.simplify(hamiltonian_action(coherent) + wd^(-grade) * dissipative)::Liouvillian
end

function native_validation_matrix_equal(left, right)
  size(left) == size(right) || return false
  return all(iszero(SQA.simplify(left[index] - right[index])) for index in eachindex(left))
end

function native_validation_seed_amplitudes(native)
  return [amplitude.seed.amplitude for amplitude in native.amplitudes]
end

function native_validation_prefix_equal(low, high, retained_grades)
  for grade in retained_grades
    difference = low.coefficients[grade + 1] - high.coefficients[grade + 1]
    all(iszero(SQA.simplify(difference[harmonic])) for harmonic in keys(difference)) ||
      return false
  end
  return true
end

@testset "native GKSL expansion matches Liouvillian van Vleck on the driven qubit" begin
  pauli = PauliSpace(:native_gksl_validation_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  bright = σy + σz
  @variables ω_native_qubit::Real t_native_qubit::Real Ω_native_qubit::Real

  ω = ω_native_qubit
  t = t_native_qubit
  H = Ω_native_qubit * cos(ω * t) * σx
  physical_channels = (collapse(bright),)

  native_order2 = FE.native_gksl_expansion(H, ω, t, 2, physical_channels)
  native = FE.native_gksl_expansion(H, ω, t, 3, physical_channels)
  raw = floquet_expansion(H, ω, t, VanVleck(), 3; channels=physical_channels)

  @test native_validation_prefix_equal(
    only(native_order2.amplitudes), only(native.amplitudes), 0:1
  )

  # This is the driven-qubit fixture of the Gram and spectral completion tests. The physical
  # collapse amplitude is static, so R_{m!=0}=0 and the second-order static similarity
  # [B_R^(2), H_0] vanishes. Amplitude transport and Liouvillian van Vleck therefore agree
  # coefficient by coefficient throughout the one-dissipator sector.
  for grade in 0:2
    @test native_validation_liouvillian_zero(
      effective_component(raw, grade) -
      native_validation_retained_component(native, ω, grade),
    )
  end
end

@testset "native GKSL expansion is covariant on the full-rank two-channel fixture" begin
  fock = FockSpace(:native_gksl_validation_full_rank)
  a = Destroy(fock, :a)
  frame = DissipativeFrame(a, a^2)
  @variables ω_native_full::Real t_native_full::Real

  ω = ω_native_full
  t = t_native_full
  H = 0 * a
  first = collapse(a + a^2)
  second = collapse(a + im * a^2)
  physical_channels = (first, second)

  native = FE.native_gksl_expansion(H, ω, t, 1, physical_channels)
  raw = floquet_expansion(H, ω, t, VanVleck(), 1; channels=physical_channels)
  gram = positive_completion(raw, Gram(), frame)

  @test native_validation_liouvillian_zero(native.generator - effective_generator(raw))
  @test native_validation_matrix_equal(
    kossakowski(native.generator, frame), kossakowski(gram)
  )

  # Reordering microscopic channels permutes the provenance but not the physical Kraus sum.
  reversed = FE.native_gksl_expansion(H, ω, t, 1, (second, first))
  @test native_validation_liouvillian_zero(native.generator - reversed.generator)
  @test [amplitude.seed.reference.index for amplitude in native.amplitudes] == [1, 2]
  @test native_validation_seed_amplitudes(reversed) ==
    reverse(native_validation_seed_amplitudes(native))
end

@testset "native GKSL expansion agrees with Gram completion across a frame congruence" begin
  fock = FockSpace(:native_gksl_validation_congruence)
  a = Destroy(fock, :a)
  native_frame = DissipativeFrame(a, a^2)
  g1 = a + (1 // 2) * a^2
  g2 = 2a - a^2
  transformed_frame = DissipativeFrame(g1, g2)
  channel = 2g1 + 3g2
  @variables ω_native_congruence::Real t_native_congruence::Real

  ω = ω_native_congruence
  t = t_native_congruence
  H = 0 * a
  physical_channels = (collapse(channel),)

  native = FE.native_gksl_expansion(H, ω, t, 1, physical_channels)
  raw = floquet_expansion(H, ω, t, VanVleck(), 1; channels=physical_channels)
  gram_native = positive_completion(raw, Gram(), native_frame)
  gram_transformed = positive_completion(raw, Gram(), transformed_frame)

  @test native_validation_liouvillian_zero(native.generator - effective_generator(raw))
  @test native_validation_matrix_equal(
    kossakowski(native.generator, native_frame), kossakowski(gram_native)
  )
  @test native_validation_matrix_equal(
    kossakowski(native.generator, transformed_frame), kossakowski(gram_transformed)
  )
end

@testset "native GKSL expansion agrees with completions on Kerr number-selective loss" begin
  fock = FockSpace(:native_gksl_validation_number_selective)
  a = Destroy(fock, :a)
  number_selective = a' * a^2
  frame = DissipativeFrame(number_selective)
  @variables ω_native_ns::Real t_native_ns::Real K_native_ns::Real γ_native_ns::Real

  ω = ω_native_ns
  t = t_native_ns
  H = K_native_ns * a'^2 * a^2
  physical_channels = (jump(number_selective, γ_native_ns),)

  native = FE.native_gksl_expansion(H, ω, t, 1, physical_channels)
  raw = floquet_expansion(H, ω, t, VanVleck(), 1; channels=physical_channels)
  gram = positive_completion(raw, Gram(), frame)
  spectral = positive_completion(raw, Spectral(), frame)

  @test native_validation_liouvillian_zero(native.generator - effective_generator(raw))
  native_matrix = kossakowski(native.generator, frame)
  @test native_validation_matrix_equal(native_matrix, kossakowski(gram))
  @test native_validation_matrix_equal(native_matrix, kossakowski(spectral))
end

@testset "native GKSL expansion turns driven two-photon loss into one-photon loss" begin
  fock = FockSpace(:native_gksl_validation_kerr)
  a = Destroy(fock, :a)
  @variables ω_native_kerr::Real t_native_kerr::Real ε_native_kerr::Real
  @variables Δ_native_kerr::Real K_native_kerr::Real

  ω = ω_native_kerr
  t = t_native_kerr
  ε = ε_native_kerr
  H = Δ_native_kerr * a' * a + K_native_kerr * a'^2 * a^2 + ε * cos(ω * t) * (a + a')

  physical_channels = (collapse(a^2),)
  native = FE.native_gksl_expansion(H, ω, t, 3, physical_channels)
  raw = floquet_expansion(H, ω, t, VanVleck(), 3; channels=physical_channels)

  # Static two-photon loss again has B_R^(2)=0, so the second-order retained coefficient
  # is a direct cutoff-free comparison with Liouvillian van Vleck.
  retained = native_validation_retained_component(native, ω, 2)
  @test native_validation_liouvillian_zero(effective_component(raw, 2) - retained)

  # The drive micromotion displaces the mode, a -> a + α(t) with |α(t)| = (ε/ω)|sin ωt|,
  # so the transported amplitude a^2 -> a^2 + 2α(t) a + α(t)^2 has one-photon harmonics
  # ∓(ε/ω) a at m = ±1. Squaring them gives the rate 4 * mean|α|^2 = 2ε^2/ω^2 without a Fock
  # cutoff.
  one_photon_rate = 2 * ε^2 / ω^2
  @test native_validation_matrix_equal(
    kossakowski(retained, DissipativeFrame(a, a^2)),
    [
      convert(SQA.CNum, one_photon_rate) convert(SQA.CNum, 0)
      convert(SQA.CNum, 0) convert(SQA.CNum, 0)
    ],
  )
end

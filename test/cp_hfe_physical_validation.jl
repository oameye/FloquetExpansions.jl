using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const FE = FloquetExpansions
const SQA = SecondQuantizedAlgebra

cp_validation_liouvillian_zero(L::Liouvillian) = iszero(FE.canonical_liouvillian(L))

function cp_validation_retained_component(reconstruction, grade::Int)
  coherent = getfield(reconstruction.coherent, :effective_components)[grade + 1]
  dissipative = FE.cp_dissipative_component(reconstruction.amplitudes, grade)
  return SQA.simplify(hamiltonian_action(coherent) + dissipative)::Liouvillian
end

function cp_validation_matrix_equal(left, right)
  size(left) == size(right) || return false
  return all(iszero(SQA.simplify(left[index] - right[index])) for index in eachindex(left))
end

function cp_validation_matrix_hermitian(matrix)
  size(matrix, 1) == size(matrix, 2) || return false
  return all(
    iszero(SQA.simplify(matrix[row, column] - conj(matrix[column, row]))) for
    row in axes(matrix, 1), column in axes(matrix, 2)
  )
end

function cp_validation_prefix_equal(low, high, retained_grades)
  for grade in retained_grades
    difference = low.coefficients[grade + 1] - high.coefficients[grade + 1]
    all(iszero(SQA.simplify(difference[harmonic])) for harmonic in keys(difference)) ||
      return false
  end
  return true
end

@testset "existing driven-qubit completion benchmark anchors native CP-HFE" begin
  pauli = PauliSpace(:cp_hfe_validation_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  bright = σy + σz
  dark = σy - σz
  cartesian = DissipativeFrame(σx, σy, σz)
  adapted = DissipativeFrame(bright, dark)
  @variables ω_cp_qubit::Real t_cp_qubit::Real Ω_cp_qubit::Real

  ω = ω_cp_qubit
  t = t_cp_qubit
  H = Ω_cp_qubit * cos(ω * t) * σx
  physical_channels = (collapse(bright),)

  native_order2 = FE.cp_hfe_reconstruction(H, ω, t, 2, physical_channels)
  native = FE.cp_hfe_reconstruction(H, ω, t, 3, physical_channels)
  raw = floquet_expansion(H, ω, t, VanVleck(), 3; channels=physical_channels)
  raw_components = getfield(raw, :effective_components)

  @test cp_validation_prefix_equal(
    only(native_order2.amplitudes), only(native.amplitudes), 0:1
  )

  # This is exactly the existing Gram/Spectral driven-qubit fixture. The physical collapse
  # amplitude is static, so R_{m!=0}=0 and the certified B_R^(2) frame correction vanishes.
  # Native amplitude transport and direct Liouvillian Van Vleck therefore agree coefficient by
  # coefficient throughout the complete one-dissipator sector.
  for grade in 0:2
    @test cp_validation_liouvillian_zero(
      raw_components[grade + 1] - cp_validation_retained_component(native, grade)
    )
  end

  gram_cartesian = positive_completion(raw, Gram(), cartesian)
  gram_adapted = positive_completion(raw, Gram(), adapted)
  spectral_adapted = positive_completion(raw, Spectral(), adapted)
  for grade in 0:2
    @test effective_component(gram_cartesian, grade) == effective_component(raw, grade)
    @test effective_component(gram_adapted, grade) == effective_component(raw, grade)
    @test effective_component(spectral_adapted, grade) == effective_component(raw, grade)
  end

  native_cartesian = kossakowski(native.generator, cartesian)
  native_adapted = kossakowski(native.generator, adapted)
  @test cp_validation_matrix_hermitian(native_cartesian)
  @test cp_validation_matrix_hermitian(native_adapted)
  @test cp_validation_matrix_hermitian(kossakowski(gram_cartesian))
  @test cp_validation_matrix_hermitian(kossakowski(gram_adapted))
  @test cp_validation_matrix_hermitian(kossakowski(spectral_adapted))
end

@testset "existing full-rank two-channel fixture preserves microscopic channel covariance" begin
  fock = FockSpace(:cp_hfe_validation_full_rank)
  a = Destroy(fock, :a)
  frame = DissipativeFrame(a, a^2)
  @variables ω_cp_full::Real t_cp_full::Real

  ω = ω_cp_full
  t = t_cp_full
  H = 0 * a
  first = collapse(a + a^2)
  second = collapse(a + im * a^2)
  physical_channels = (first, second)

  native = FE.cp_hfe_reconstruction(H, ω, t, 1, physical_channels)
  raw = floquet_expansion(H, ω, t, VanVleck(), 1; channels=physical_channels)
  gram = positive_completion(raw, Gram(), frame)

  @test cp_validation_liouvillian_zero(native.generator - effective_generator(raw))
  @test cp_validation_matrix_equal(kossakowski(native.generator, frame), kossakowski(gram))

  # Reordering microscopic channels changes provenance order but not the physical Kraus sum.
  reversed = FE.cp_hfe_reconstruction(H, ω, t, 1, (second, first))
  @test cp_validation_liouvillian_zero(native.generator - reversed.generator)
  @test native.amplitudes[1].seed.reference.index == 1
  @test native.amplitudes[2].seed.reference.index == 2
  @test reversed.amplitudes[1].seed.reference.index == 1
  @test reversed.amplitudes[2].seed.reference.index == 2
end

@testset "existing rational frame-congruence fixture stays physical in native CP-HFE" begin
  fock = FockSpace(:cp_hfe_validation_congruence)
  a = Destroy(fock, :a)
  native_frame = DissipativeFrame(a, a^2)
  g1 = a + (1 // 2) * a^2
  g2 = 2a - a^2
  transformed_frame = DissipativeFrame(g1, g2)
  channel = 2g1 + 3g2
  @variables ω_cp_congruence::Real t_cp_congruence::Real

  ω = ω_cp_congruence
  t = t_cp_congruence
  H = 0 * a
  physical_channels = (collapse(channel),)

  native = FE.cp_hfe_reconstruction(H, ω, t, 1, physical_channels)
  raw = floquet_expansion(H, ω, t, VanVleck(), 1; channels=physical_channels)
  gram_native = positive_completion(raw, Gram(), native_frame)
  gram_transformed = positive_completion(raw, Gram(), transformed_frame)

  @test cp_validation_liouvillian_zero(native.generator - effective_generator(raw))
  @test cp_validation_matrix_equal(
    kossakowski(native.generator, native_frame), kossakowski(gram_native)
  )
  @test cp_validation_matrix_equal(
    kossakowski(native.generator, transformed_frame), kossakowski(gram_transformed)
  )
  @test iszero(
    SQA.simplify(
      hamiltonian(native.generator, native_frame) -
      hamiltonian(native.generator, transformed_frame)
    )
  )
end

@testset "existing Kerr number-selective-loss fixture agrees across CP constructions" begin
  fock = FockSpace(:cp_hfe_validation_number_selective)
  a = Destroy(fock, :a)
  number_selective = a' * a^2
  frame = DissipativeFrame(number_selective)
  @variables ω_cp_ns::Real t_cp_ns::Real K_cp_ns::Real γ_cp_ns::Real

  ω = ω_cp_ns
  t = t_cp_ns
  H = K_cp_ns * a'^2 * a^2
  physical_channels = (jump(number_selective, γ_cp_ns),)

  native = FE.cp_hfe_reconstruction(H, ω, t, 1, physical_channels)
  raw = floquet_expansion(H, ω, t, VanVleck(), 1; channels=physical_channels)
  gram = positive_completion(raw, Gram(), frame)
  spectral = positive_completion(raw, Spectral(), frame)

  @test cp_validation_liouvillian_zero(native.generator - effective_generator(raw))
  native_matrix = kossakowski(native.generator, frame)
  @test cp_validation_matrix_equal(native_matrix, kossakowski(gram))
  @test cp_validation_matrix_equal(native_matrix, kossakowski(spectral))
  @test cp_validation_matrix_hermitian(native_matrix)
end

@testset "driven Kerr resonator generates one-photon loss from two-photon loss" begin
  fock = FockSpace(:cp_hfe_validation_kerr)
  a = Destroy(fock, :a)
  @variables ω_cp_kerr::Real t_cp_kerr::Real Δ_cp_kerr::Real K_cp_kerr::Real ε_cp_kerr::Real

  ω = ω_cp_kerr
  t = t_cp_kerr
  number = a' * a
  H = Δ_cp_kerr * number + K_cp_kerr * a'^2 * a^2 + ε_cp_kerr * cos(ω * t) * (a + a')

  physical_channels = (collapse(a^2),)
  native = FE.cp_hfe_reconstruction(H, ω, t, 3, physical_channels)
  raw = floquet_expansion(H, ω, t, VanVleck(), 3; channels=physical_channels)

  # Static two-photon loss again has B_R^(2)=0. The second-order retained map is therefore a
  # direct cutoff-free comparison with generic Liouvillian Van Vleck.
  @test cp_validation_liouvillian_zero(
    getfield(raw, :effective_components)[3] - cp_validation_retained_component(native, 2)
  )

  # The first coherent kick dresses L=a^2 by i[K^(1),L]. A linear coherent drive creates an `a`
  # amplitude. Squaring that transported first-order amplitude produces one-photon loss at second
  # order without introducing a Fock cutoff.
  first_dressed = only(native.amplitudes).coefficients[2]
  generated = zero(Liouvillian)
  for harmonic in keys(first_dressed)
    component = SQA.simplify(first_dressed[harmonic])
    iszero(component) && continue
    generated += dissipator(component)
  end

  one_photon = kossakowski(SQA.simplify(generated), DissipativeFrame(a))
  @test size(one_photon) == (1, 1)
  @test !iszero(SQA.simplify(one_photon[1, 1]))
end

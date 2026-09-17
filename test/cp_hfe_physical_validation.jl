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

function cp_validation_matrix_hermitian(matrix)
  size(matrix, 1) == size(matrix, 2) || return false
  return all(
    iszero(SQA.simplify(matrix[row, column] - conj(matrix[column, row]))) for
    row in axes(matrix, 1), column in axes(matrix, 2)
  )
end

@testset "driven qubit validates the physical one-dissipator sector" begin
  pauli = PauliSpace(:cp_hfe_validation_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  bright = σy + σz
  dark = σy - σz
  frame = DissipativeFrame(bright, dark)
  @variables ω_cp_qubit::Real t_cp_qubit::Real Δ_cp_qubit::Real Ω_cp_qubit::Real

  ω = ω_cp_qubit
  t = t_cp_qubit
  H = Δ_cp_qubit * σz + Ω_cp_qubit * cos(ω * t) * σx
  physical_channels = (collapse(bright),)

  native_order2 = FE.cp_hfe_reconstruction(H, ω, t, 2, physical_channels)
  native_order3 = FE.cp_hfe_reconstruction(H, ω, t, 3, physical_channels)

  # Requesting a higher truncation order must not alter already-retained physical amplitudes.
  low_series = only(native_order2.amplitudes)
  high_series = only(native_order3.amplitudes)
  for grade in 0:1
    difference = low_series.coefficients[grade + 1] - high_series.coefficients[grade + 1]
    @test all(iszero(SQA.simplify(difference[harmonic])) for harmonic in keys(difference))
  end

  raw = floquet_expansion(H, ω, t, VanVleck(), 3; channels=physical_channels)
  raw_components = getfield(raw, :effective_components)

  # Static physical loss carries no nonzero dissipative Fourier harmonic, hence B_R^(2)=0.
  # The complete one-dissipator sector therefore agrees directly in this physical frame.
  for grade in 0:2
    native_component = cp_validation_retained_component(native_order3, grade)
    @test cp_validation_liouvillian_zero(raw_components[grade + 1] - native_component)
  end

  gram = positive_completion(raw, Gram(), frame)
  spectral = positive_completion(raw, Spectral(), frame)
  for grade in 0:2
    @test effective_component(gram, grade) == effective_component(raw, grade)
    @test effective_component(spectral, grade) == effective_component(raw, grade)
  end

  native_kossakowski = kossakowski(native_order3.generator, frame)
  @test cp_validation_matrix_hermitian(native_kossakowski)
  @test cp_validation_matrix_hermitian(kossakowski(gram))
  @test cp_validation_matrix_hermitian(kossakowski(spectral))
end

@testset "driven Kerr resonator generates one-photon loss from two-photon loss" begin
  fock = FockSpace(:cp_hfe_validation_kerr)
  a = Destroy(fock, :a)
  @variables ω_cp_kerr::Real t_cp_kerr::Real Δ_cp_kerr::Real K_cp_kerr::Real ε_cp_kerr::Real

  ω = ω_cp_kerr
  t = t_cp_kerr
  number = a' * a
  H =
    Δ_cp_kerr * number +
    K_cp_kerr * a'^2 * a^2 +
    ε_cp_kerr * cos(ω * t) * (a + a')

  physical_channels = (collapse(a^2),)
  native = FE.cp_hfe_reconstruction(H, ω, t, 3, physical_channels)
  raw = floquet_expansion(H, ω, t, VanVleck(), 3; channels=physical_channels)

  # With a static physical dissipator the retained one-dissipator sector must again agree with
  # direct Liouvillian Van Vleck through second order, now entirely in symbolic bosonic algebra.
  @test cp_validation_liouvillian_zero(
    getfield(raw, :effective_components)[3] - cp_validation_retained_component(native, 2)
  )

  # The first coherent kick dresses L=a^2 by i[K^(1),L]. For a linear coherent drive this
  # creates an `a` direction. Squaring that transported first-order amplitude therefore produces
  # a genuine one-photon dissipator at second order, without introducing a Fock cutoff.
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

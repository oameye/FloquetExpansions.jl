using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA_GKLS_VALIDATION = SecondQuantizedAlgebra

function gkls_validation_zero(value)
  return iszero(SQA_GKLS_VALIDATION.simplify(value))
end

function gkls_validation_matrix_equal(left, right)
  size(left) == size(right) || return false
  return all(
    iszero(SQA_GKLS_VALIDATION.simplify(left[index] - right[index])) for
    index in eachindex(left)
  )
end

function gkls_validation_matrix_hermitian(matrix)
  size(matrix, 1) == size(matrix, 2) || return false
  return all(
    iszero(SQA_GKLS_VALIDATION.simplify(matrix[row, column] - conj(matrix[column, row])))
    for row in axes(matrix, 1), column in axes(matrix, 2)
  )
end

function gkls_validation_pair(H, ω, t, order; channels)
  hd = floquet_expansion(H, ω, t, VanVleck(; algorithm=HoriDeprit()), order; channels)
  bf = floquet_expansion(H, ω, t, VanVleck(; algorithm=BlochFeshbach()), order; channels)
  hd_raw = floquet_expansion(
    H,
    ω,
    t,
    VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))),
    order;
    channels,
  )
  bf_raw = floquet_expansion(
    H,
    ω,
    t,
    VanVleck(; algorithm=BlochFeshbach(; complete_positive=Val(false))),
    order;
    channels,
  )
  return (; hd, bf, hd_raw, bf_raw)
end

function gkls_validation_certify(pair, order)
  @test factorization(pair.hd) isa FloquetExpansions.GramFactorization
  @test factorization(pair.bf) isa FloquetExpansions.GramFactorization

  hd_expected = positive_completion(pair.hd_raw, Gram())
  bf_expected = positive_completion(pair.bf_raw, Gram())
  @test gkls_validation_zero(
    effective_generator(pair.hd) - effective_generator(hd_expected)
  )
  @test gkls_validation_zero(
    effective_generator(pair.bf) - effective_generator(bf_expected)
  )
  @test gkls_validation_zero(effective_generator(pair.hd) - effective_generator(pair.bf))

  @test gkls_validation_zero(
    liouvillian(hamiltonian(pair.hd); channels=channels(pair.hd)) -
    effective_generator(pair.hd),
  )
  @test gkls_validation_zero(
    liouvillian(hamiltonian(pair.bf); channels=channels(pair.bf)) -
    effective_generator(pair.bf),
  )

  for grade in 0:(order - 1)
    @test gkls_validation_zero(
      effective_component(pair.hd_raw, grade) - effective_component(pair.bf_raw, grade)
    )
    @test gkls_validation_zero(
      effective_component(pair.hd, grade) - effective_component(pair.hd_raw, grade)
    )
    @test gkls_validation_zero(
      effective_component(pair.bf, grade) - effective_component(pair.bf_raw, grade)
    )
  end
  for grade in 1:(order - 1)
    @test gkls_validation_zero(
      micromotion(pair.hd_raw, grade) - micromotion(pair.bf_raw, grade)
    )
    @test gkls_validation_zero(
      micromotion(pair.hd, grade) - micromotion(pair.hd_raw, grade)
    )
    @test gkls_validation_zero(
      micromotion(pair.bf, grade) - micromotion(pair.bf_raw, grade)
    )
  end
  return nothing
end

@testset "GKLS HFE physical validation: driven qubit" begin
  pauli = PauliSpace(:gkls_hfe_validation_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  bright = σy + σz
  dark = σy - σz
  cartesian = DissipativeFrame(σx, σy, σz)
  adapted = DissipativeFrame(bright, dark)
  @variables ω_gkls_qubit::Real t_gkls_qubit::Real Ω_gkls_qubit::Real

  H = Ω_gkls_qubit * cos(ω_gkls_qubit * t_gkls_qubit) * σx
  physical_channels = (collapse(bright),)
  pair = gkls_validation_pair(H, ω_gkls_qubit, t_gkls_qubit, 3; channels=physical_channels)
  gkls_validation_certify(pair, 3)

  gram_cartesian = positive_completion(pair.hd_raw, Gram(), cartesian)
  gram_adapted = positive_completion(pair.hd_raw, Gram(), adapted)
  spectral_adapted = positive_completion(pair.hd_raw, Spectral(), adapted)
  for grade in 0:2
    @test gkls_validation_zero(
      effective_component(gram_cartesian, grade) - effective_component(pair.hd_raw, grade)
    )
    @test gkls_validation_zero(
      effective_component(gram_adapted, grade) - effective_component(pair.hd_raw, grade)
    )
    @test gkls_validation_zero(
      effective_component(spectral_adapted, grade) - effective_component(pair.hd_raw, grade)
    )
  end
  @test gkls_validation_matrix_hermitian(
    kossakowski(effective_generator(pair.hd), cartesian)
  )
  @test gkls_validation_matrix_hermitian(kossakowski(effective_generator(pair.hd), adapted))
end

@testset "GKLS HFE physical validation: full-rank two-channel bosonic model" begin
  fock = FockSpace(:gkls_hfe_validation_full_rank)
  a = Destroy(fock, :a)
  frame = DissipativeFrame(a, a^2)
  @variables ω_gkls_full::Real t_gkls_full::Real

  physical_channels = (collapse(a + a^2), collapse(a + im * a^2))
  pair = gkls_validation_pair(
    0 * a, ω_gkls_full, t_gkls_full, 1; channels=physical_channels
  )
  gkls_validation_certify(pair, 1)

  explicit_gram = positive_completion(pair.hd_raw, Gram(), frame)
  @test factorization(explicit_gram) isa FloquetExpansions.GramFactorization
  @test gkls_validation_zero(
    effective_component(explicit_gram, 0) - effective_component(pair.hd_raw, 0)
  )
  @test gkls_validation_matrix_equal(
    kossakowski_component(explicit_gram, 0), kossakowski_component(pair.hd_raw, frame, 0)
  )
  @test gkls_validation_zero(
    liouvillian(hamiltonian(explicit_gram); channels=channels(explicit_gram)) -
    effective_generator(explicit_gram),
  )
  @test gkls_validation_matrix_hermitian(kossakowski(explicit_gram))

  # At order one the Gram factor contains only the retained leading form, so automatic and fixed
  # frames represent the same completed generator. Compare in common GKLS coordinates: direct
  # subtraction of the two nonunitarily related symbolic Liouvillian representations is not a
  # reliable canonicalization oracle for SQA.
  @test gkls_validation_zero(hamiltonian(pair.hd) - hamiltonian(explicit_gram))
  @test gkls_validation_matrix_equal(
    kossakowski(effective_generator(pair.hd), frame), kossakowski(explicit_gram)
  )
end

@testset "GKLS HFE physical validation: rational frame congruence" begin
  fock = FockSpace(:gkls_hfe_validation_congruence)
  a = Destroy(fock, :a)
  native_frame = DissipativeFrame(a, a^2)
  g1 = a + (1 // 2) * a^2
  g2 = 2a - a^2
  transformed_frame = DissipativeFrame(g1, g2)
  channel = 2g1 + 3g2
  @variables ω_gkls_congruence::Real t_gkls_congruence::Real

  pair = gkls_validation_pair(
    0 * a, ω_gkls_congruence, t_gkls_congruence, 1; channels=(collapse(channel),)
  )
  gkls_validation_certify(pair, 1)

  native_matrix = kossakowski(effective_generator(pair.hd), native_frame)
  transformed_matrix = kossakowski(effective_generator(pair.hd), transformed_frame)
  @test gkls_validation_matrix_hermitian(native_matrix)
  @test gkls_validation_matrix_hermitian(transformed_matrix)
  @test gkls_validation_zero(
    hamiltonian(effective_generator(pair.hd), native_frame) -
    hamiltonian(effective_generator(pair.hd), transformed_frame),
  )
end

@testset "GKLS HFE physical validation: Kerr number-selective loss" begin
  fock = FockSpace(:gkls_hfe_validation_number_selective)
  a = Destroy(fock, :a)
  number_selective = a' * a^2
  frame = DissipativeFrame(number_selective)
  @variables ω_gkls_ns::Real t_gkls_ns::Real K_gkls_ns::Real γ_gkls_ns::Real

  H = K_gkls_ns * a'^2 * a^2
  pair = gkls_validation_pair(
    H, ω_gkls_ns, t_gkls_ns, 1; channels=(jump(number_selective, γ_gkls_ns),)
  )
  gkls_validation_certify(pair, 1)

  gram = positive_completion(pair.hd_raw, Gram(), frame)
  spectral = positive_completion(pair.hd_raw, Spectral(), frame)
  @test gkls_validation_matrix_equal(
    kossakowski(effective_generator(pair.hd), frame), kossakowski(gram)
  )
  @test gkls_validation_matrix_equal(kossakowski(gram), kossakowski(spectral))
end

@testset "GKLS HFE physical validation: driven Kerr with two-photon loss" begin
  fock = FockSpace(:gkls_hfe_validation_kerr)
  a = Destroy(fock, :a)
  frame = DissipativeFrame(a, a^2)
  @variables ω_gkls_kerr::Real t_gkls_kerr::Real Δ_gkls_kerr::Real K_gkls_kerr::Real ε_gkls_kerr::Real

  number = a' * a
  H =
    Δ_gkls_kerr * number +
    K_gkls_kerr * a'^2 * a^2 +
    ε_gkls_kerr * cos(ω_gkls_kerr * t_gkls_kerr) * (a + a')
  pair = gkls_validation_pair(H, ω_gkls_kerr, t_gkls_kerr, 3; channels=(collapse(a^2),))
  gkls_validation_certify(pair, 3)

  kossakowski_matrix = kossakowski(effective_generator(pair.hd), frame)
  @test gkls_validation_matrix_hermitian(kossakowski_matrix)

  d0 = kossakowski_component(pair.hd_raw, frame, 0)
  d1 = kossakowski_component(pair.hd_raw, frame, 1)
  d2 = kossakowski_component(pair.hd_raw, frame, 2)
  @test gkls_validation_zero(d0[1, 1])
  @test gkls_validation_zero(d1[1, 1])
  @test !iszero(SQA_GKLS_VALIDATION.simplify(d2[1, 1]))
end

using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

function analytic_matrix(rows::Tuple...)
  isempty(rows) && return Matrix{SQA.CNum}(undef, 0, 0)
  columns = length(first(rows))
  all(length(row) == columns for row in rows) || throw(DimensionMismatch())
  result = fill(convert(SQA.CNum, 0), length(rows), columns)
  for row_index in eachindex(rows), column_index in 1:columns
    result[row_index, column_index] = convert(SQA.CNum, rows[row_index][column_index])
  end
  return result
end

function analytic_matrix_equal(left, right)
  size(left) == size(right) || return false
  return all(iszero(SQA.simplify(left[index] - right[index])) for index in eachindex(left))
end

@testset "analytical driven-qubit root and HCM completions" begin
  pauli = PauliSpace(:cp_analytic_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  σminus = (1 // 2) * (σx - im * σy)
  σplus = (1 // 2) * (σx + im * σy)

  cartesian = DissipativeFrame(σx, σy, σz)
  spectral_frame = DissipativeFrame(σminus, σplus, σz)

  @variables ω::Real t::Real E::Real γ::Real
  H = (1 // 2) * σz + E * cos(ω * t) * σx
  expansion = floquet_expansion(H, ω, t, VanVleck(), 3; channels=(jump(σminus, γ),))

  z2 = E^2 / ω^2
  z4 = z2^2

  expected_cartesian = analytic_matrix(
    (γ / 4, (im * γ / 4) * (1 - z2), 0),
    ((-im * γ / 4) * (1 - z2), (γ / 4) * (1 - 2z2), 0),
    (0, 0, (γ / 2) * z2),
  )
  @test analytic_matrix_equal(kossakowski(expansion, cartesian), expected_cartesian)

  # This is the same analytical (u,v,e_z) basis as in the derivation, rescaled by
  # u -> σ₋ = u/sqrt(2), v -> σ₊ = v/sqrt(2). The rational rescaling avoids
  # introducing symbolic radicals into the public dissipative-frame extraction.
  expected_spectral_target = analytic_matrix(
    (γ * (1 - z2), (γ / 2) * z2, 0), ((γ / 2) * z2, 0, 0), (0, 0, (γ / 2) * z2)
  )
  @test analytic_matrix_equal(
    kossakowski(expansion, spectral_frame), expected_spectral_target
  )

  @test iszero(SQA.simplify(hamiltonian_component(expansion, 0) - (1 // 2) * σz))
  @test iszero(hamiltonian_component(expansion, 1))
  @test iszero(SQA.simplify(hamiltonian_component(expansion, 2) + (z2 / 2) * σz))

  gram = positive_completion(expansion, Gram(), spectral_frame)
  spectral = positive_completion(expansion, Spectral(), spectral_frame)

  expected_gram = analytic_matrix(
    (γ * (1 - z2 + z4 / 4), (γ / 2) * z2 - (γ / 4) * z4, 0),
    ((γ / 2) * z2 - (γ / 4) * z4, (γ / 4) * z4, 0),
    (0, 0, (γ / 2) * z2),
  )
  @test analytic_matrix_equal(kossakowski(gram), expected_gram)

  hcm_active = 1 - z2 + z4 / 4
  hcm_denominator = 1 + z4 / 4
  hcm_prefactor = γ / hcm_denominator
  expected_hcm = analytic_matrix(
    (hcm_prefactor * hcm_active, hcm_prefactor * (z2 / 2) * hcm_active, 0),
    (hcm_prefactor * (z2 / 2) * hcm_active, hcm_prefactor * (z4 / 4) * hcm_active, 0),
    (0, 0, (γ / 2) * z2),
  )
  @test analytic_matrix_equal(kossakowski(spectral), expected_hcm)
  @test !analytic_matrix_equal(kossakowski(gram), kossakowski(spectral))

  for order in 0:2
    @test effective_component(gram, order) == effective_component(expansion, order)
    @test effective_component(spectral, order) == effective_component(expansion, order)
    @test analytic_matrix_equal(
      kossakowski_component(gram, order),
      kossakowski_component(expansion, spectral_frame, order),
    )
    @test analytic_matrix_equal(
      kossakowski_component(spectral, order),
      kossakowski_component(expansion, spectral_frame, order),
    )
  end
  @test micromotion(gram) == micromotion(expansion)
  @test micromotion(spectral) == micromotion(expansion)

  gram_data = factorization(gram)
  spectral_data = factorization(spectral)
  @test gram_data.onsets == [0, 1]
  @test spectral_data.onsets == [0, -1, 2]
  @test spectral_data.puiseux == [false, false, false]

  comparison_prefactor = γ * z4 / (4 * hcm_denominator)
  expected_difference = analytic_matrix(
    (comparison_prefactor * hcm_active, comparison_prefactor * (1 - z4 / 4), 0),
    (comparison_prefactor * (1 - z4 / 4), comparison_prefactor * z2, 0),
    (0, 0, 0),
  )
  @test analytic_matrix_equal(
    kossakowski(gram) - kossakowski(spectral), expected_difference
  )
end

@testset "analytical modulated one- and two-photon loss closure" begin
  fock = FockSpace(:cp_analytic_bosonic_loss)
  a = Destroy(fock, :a)
  na = a' * a^2
  frame = DissipativeFrame(a, a^2, na)
  @variables ω::Real κ1::Real κ2::Real r1::Real r2::Real

  u = κ1 * r1 / 2
  v = κ2 * r2 / 2
  χ = κ1 * κ2 * r1 * r2 / 2
  D1 = dissipator(a)
  D2 = dissipator(a^2)

  # Quarter-cycle reservoir phase, φ = π/2, in the analytical model
  # κ1[1+r1 cos(ωt)]D[a] + κ2[1+r2 cos(ωt-φ)]D[a²].
  generator = PeriodicGenerator(
    Dict(0 => κ1 * D1 + κ2 * D2, 1 => u * D1 + im * v * D2, -1 => u * D1 - im * v * D2), ω
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)
  c = χ / ω

  expected_leading = analytic_matrix((κ1, 0, 0), (0, κ2, 0), (0, 0, 0))
  expected_first_order = analytic_matrix((0, 0, c), (0, -2c, 0), (c, 0, 0))
  @test analytic_matrix_equal(kossakowski_component(expansion, frame, 0), expected_leading)
  @test analytic_matrix_equal(
    kossakowski_component(expansion, frame, 1), expected_first_order
  )

  expected_direct = analytic_matrix((κ1, 0, c), (0, κ2 - 2c, 0), (c, 0, 0))
  @test analytic_matrix_equal(kossakowski(expansion, frame), expected_direct)

  gram = positive_completion(expansion, Gram(), frame)
  expected_completion = analytic_matrix(
    (κ1, 0, c), (0, κ2 - 2c + c^2 / κ2, 0), (c, 0, c^2 / κ1)
  )
  @test analytic_matrix_equal(kossakowski(gram), expected_completion)

  closure = kossakowski(gram) - kossakowski(expansion, frame)
  expected_closure = analytic_matrix((0, 0, 0), (0, c^2 / κ2, 0), (0, 0, c^2 / κ1))
  @test analytic_matrix_equal(closure, expected_closure)

  gram_data = factorization(gram)
  @test gram_data.onsets == [0, 0]
  @test length(gram_data.amplitudes) == 2
  @test liouvillian(hamiltonian(gram); channels=channels(gram)) == effective_generator(gram)
  @test micromotion(gram) == micromotion(expansion)
end

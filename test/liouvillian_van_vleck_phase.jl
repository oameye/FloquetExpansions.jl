using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

function phase_matrix(rows::Tuple...)
  columns = length(first(rows))
  result = fill(convert(SQA.CNum, 0), length(rows), columns)
  for row in eachindex(rows), column in 1:columns
    result[row, column] = convert(SQA.CNum, rows[row][column])
  end
  return result
end

function phase_matrix_equal(left, right)
  size(left) == size(right) || return false
  for index in eachindex(left)
    iszero(SQA.simplify(left[index] - right[index])) || return false
  end
  return true
end

@testset "Liouvillian Van Vleck phase reproduces the analytical driven qubit" begin
  pauli = PauliSpace(:liouvillian_vv_phase_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  σminus = (1 // 2) * (σx - im * σy)
  frame = DissipativeFrame(σx, σy, σz)
  @variables ω::Real t::Real E::Real γ::Real

  H = (1 // 2) * σz + E * cos(ω * t) * σx
  parsed = floquet_expansion(H, ω, t, VanVleck(), 3; channels=(jump(σminus, γ),))

  L0 = liouvillian((1 // 2) * σz; channels=(jump(σminus, γ),))
  L1 = hamiltonian_action(((1 // 2) * E) * σx)
  explicit = floquet_expansion(
    PeriodicGenerator(Dict(0 => L0, 1 => L1, -1 => L1), ω), VanVleck(), 3
  )

  for order in 0:2
    @test effective_component(parsed, order) == effective_component(explicit, order)
  end
  @test micromotion(parsed) == micromotion(explicit)

  @test iszero(hamiltonian_component(parsed, 1))
  @test iszero(
    SQA.simplify(hamiltonian_component(parsed, 2) + (E^2 / (2 * ω^2)) * σz)
  )

  expected_second = phase_matrix(
    (0, -im * γ * E^2 / (4 * ω^2), 0),
    (im * γ * E^2 / (4 * ω^2), -γ * E^2 / (2 * ω^2), 0),
    (0, 0, γ * E^2 / (2 * ω^2)),
  )
  @test phase_matrix_equal(kossakowski_component(parsed, frame, 2), expected_second)

  z2 = E^2 / ω^2
  expected_finite = phase_matrix(
    (γ / 4, (im * γ / 4) * (1 - z2), 0),
    ((-im * γ / 4) * (1 - z2), (γ / 4) * (1 - 2z2), 0),
    (0, 0, (γ / 2) * z2),
  )
  @test phase_matrix_equal(kossakowski(parsed, frame), expected_finite)
end

@testset "Liouvillian phase gives the analytical modulated-loss commutator" begin
  fock = FockSpace(:liouvillian_vv_phase_boson)
  a = Destroy(fock, :a)
  na = a' * a^2
  frame = DissipativeFrame(a, a^2, na)
  @variables ω::Real t::Real κ1::Real κ2::Real r1::Real r2::Real

  D1 = dissipator(a)
  D2 = dissipator(a^2)
  u = κ1 * r1 / 2
  v = κ2 * r2 / 2
  χ = κ1 * κ2 * r1 * r2 / 2

  one_photon_loss = κ1 * (1 + r1 * cos(ω * t)) * D1
  two_photon_loss = κ2 * (1 + r2 * sin(ω * t)) * D2
  time_domain = one_photon_loss + two_photon_loss
  parsed = floquet_expansion(time_domain, ω, t, VanVleck(), 2)
  explicit = floquet_expansion(
    PeriodicGenerator(
      Dict(
        0 => κ1 * D1 + κ2 * D2,
        1 => u * D1 + im * v * D2,
        -1 => u * D1 - im * v * D2,
      ),
      ω,
    ),
    VanVleck(),
    2,
  )

  for order in 0:1
    @test effective_component(parsed, order) == effective_component(explicit, order)
  end
  @test micromotion(parsed) == micromotion(explicit)

  expected_first = phase_matrix(
    (0, 0, χ / ω),
    (0, -2 * χ / ω, 0),
    (χ / ω, 0, 0),
  )
  @test phase_matrix_equal(kossakowski_component(parsed, frame, 1), expected_first)
end

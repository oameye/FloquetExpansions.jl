using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

symbolic_zero(expression) = iszero(SQA.simplify(expression))
same_map(A, B) = iszero(SQA.expand_completeness(A - B))
function same_periodic_map(A::PeriodicGenerator, B::PeriodicGenerator)
  return all(same_map(A[l], B[l]) for l in union(keys(A), keys(B)))
end

@testset "Liouvillian Van Vleck phase reproduces the analytical driven qubit" begin
  pauli = PauliSpace(:liouvillian_vv_phase_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  σminus = (1 // 2) * (σx - im * σy)
  frame = JumpBasis(σx, σy, σz)
  @variables ω::Real t::Real E::Real γ::Real

  H = (1 // 2) * σz + E * cos(ω * t) * σx
  parsed = floquet_expansion(
    hamiltonian_action(H) + γ * dissipator(σminus), ω, t, VanVleck(), 3
  )

  L0 = SQA.SuperOperator(hamiltonian_action((1 // 2) * σz) + γ * dissipator(σminus))
  L1 = SQA.SuperOperator(hamiltonian_action(((1 // 2) * E) * σx))
  explicit = floquet_expansion(
    PeriodicGenerator(Dict(0 => L0, 1 => L1, -1 => L1), ω), VanVleck(), 3
  )

  for order in 0:2
    @test same_map(effective_component(parsed, order), effective_component(explicit, order))
  end
  @test same_periodic_map(micromotion(parsed), micromotion(explicit))

  @test iszero(hamiltonian_component(parsed, 1))
  @test symbolic_zero(hamiltonian_component(parsed, 2) + (E^2 / (2 * ω^2)) * σz)

  d2 = kossakowski_component(parsed, frame, 2)
  @test symbolic_zero(d2[1, 1])
  @test symbolic_zero(d2[1, 2] + im * γ * E^2 / (4 * ω^2))
  @test symbolic_zero(d2[1, 3])
  @test symbolic_zero(d2[2, 1] - im * γ * E^2 / (4 * ω^2))
  @test symbolic_zero(d2[2, 2] + γ * E^2 / (2 * ω^2))
  @test symbolic_zero(d2[2, 3])
  @test symbolic_zero(d2[3, 1])
  @test symbolic_zero(d2[3, 2])
  @test symbolic_zero(d2[3, 3] - γ * E^2 / (2 * ω^2))

  z2 = E^2 / ω^2
  d = kossakowski(parsed, frame)
  @test symbolic_zero(d[1, 1] - γ / 4)
  @test symbolic_zero(d[1, 2] - (im * γ / 4) * (1 - z2))
  @test symbolic_zero(d[1, 3])
  @test symbolic_zero(d[2, 1] + (im * γ / 4) * (1 - z2))
  @test symbolic_zero(d[2, 2] - (γ / 4) * (1 - 2z2))
  @test symbolic_zero(d[2, 3])
  @test symbolic_zero(d[3, 1])
  @test symbolic_zero(d[3, 2])
  @test symbolic_zero(d[3, 3] - (γ / 2) * z2)
end

@testset "Liouvillian phase gives the analytical modulated-loss commutator" begin
  fock = FockSpace(:liouvillian_vv_phase_boson)
  a = Destroy(fock, :a)
  na = a' * a^2
  frame = JumpBasis(a, a^2, na)
  @variables ω::Real t::Real κ1::Real κ2::Real r1::Real r2::Real

  D1 = SQA.SuperOperator(dissipator(a))
  D2 = SQA.SuperOperator(dissipator(a^2))
  u = κ1 * r1 / 2
  v = κ2 * r2 / 2
  χ = κ1 * κ2 * r1 * r2 / 2

  one_photon_loss = κ1 * (1 + r1 * cos(ω * t)) * D1
  two_photon_loss = κ2 * (1 + r2 * sin(ω * t)) * D2
  time_domain = one_photon_loss + two_photon_loss
  parsed = floquet_expansion(time_domain, ω, t, VanVleck(), 2)
  explicit = floquet_expansion(
    PeriodicGenerator(
      Dict(0 => κ1 * D1 + κ2 * D2, 1 => u * D1 + im * v * D2, -1 => u * D1 - im * v * D2), ω
    ),
    VanVleck(),
    2,
  )

  for order in 0:1
    @test same_map(effective_component(parsed, order), effective_component(explicit, order))
  end
  @test same_periodic_map(micromotion(parsed), micromotion(explicit))

  d1 = kossakowski_component(parsed, frame, 1)
  @test symbolic_zero(d1[1, 1])
  @test symbolic_zero(d1[1, 2])
  @test symbolic_zero(d1[1, 3] - χ / ω)
  @test symbolic_zero(d1[2, 1])
  @test symbolic_zero(d1[2, 2] + 2χ / ω)
  @test symbolic_zero(d1[2, 3])
  @test symbolic_zero(d1[3, 1] - χ / ω)
  @test symbolic_zero(d1[3, 2])
  @test symbolic_zero(d1[3, 3])
end

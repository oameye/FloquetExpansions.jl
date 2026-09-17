using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

@testset "physical periodic jump exposes one-R frame correction and RR boundary" begin
  pauli = PauliSpace(:cp_hfe_periodic_jump)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  σminus = (σx - im * σy) / 2

  @variables ω_cp_jump::Real t_cp_jump::Real Ω_cp_jump::Real Δ_cp_jump::Real
  @variables ε_cp_jump::Real γ_cp_jump::Real
  ω = ω_cp_jump
  t = t_cp_jump
  Ω = Ω_cp_jump
  Δ = Δ_cp_jump
  ε = ε_cp_jump
  γ = γ_cp_jump

  sideband = cos(ω * t) - im * sin(ω * t)
  L = σminus + ε * sideband * σz
  H = Δ * σz + Ω * cos(ω * t) * σx
  channels = (jump(L, γ),)

  R = harmonics(liouvillian(0 * σx; channels=channels), ω, t)
  H_harmonics = harmonics(H, ω, t)
  H_map = PeriodicGenerator(
    Dict(harmonic => hamiltonian_action(H_harmonics[harmonic]) for harmonic in keys(H_harmonics)),
    ω,
  )

  B_R = zero(R[0])
  for harmonic in keys(H_map)
    iszero(harmonic) && continue
    iszero(R[-harmonic]) && continue
    B_R +=
      (1 // (2 * harmonic^2)) * SQA.commutator(H_map[harmonic], R[-harmonic])
  end
  @test !iszero(SQA.simplify(B_R))

  R_osc = PeriodicGenerator(
    Dict(harmonic => R[harmonic] for harmonic in keys(R) if !iszero(harmonic)),
    ω,
  )
  direct = floquet_expansion(R_osc, VanVleck(), 2)
  RR = zero(R[0])
  for harmonic in keys(R_osc)
    RR += (im // harmonic) * compose(R_osc[-harmonic], R_osc[harmonic])
  end

  @test !iszero(SQA.simplify(RR))
  @test iszero(SQA.simplify(ω * effective_component(direct, 1) - RR))

  native = FloquetExpansions.cp_hfe_reconstruction(0 * σx, ω, t, 2, channels)
  @test iszero(native.amplitudes[1].coefficients[2])
end

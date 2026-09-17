using Test
using FloquetExpansions
using JET: JET
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQAJet = SecondQuantizedAlgebra

@testset "JET" begin
  JET.test_package(FloquetExpansions; target_modules=(FloquetExpansions,))
end

@testset "completion optimizer stability" begin
  fock = FockSpace(:jet_completion_fock)
  a = Destroy(fock, :a)
  @variables ω::Real t::Real
  gram_frame = DissipativeFrame(a, a^2)
  gram_generator = liouvillian(0 * a; channels=(collapse(a + a^2), collapse(a + im * a^2)))
  gram_expansion = floquet_expansion(gram_generator, ω, t, VanVleck(), 1)

  JET.@test_opt target_modules=(FloquetExpansions,) positive_completion(
    gram_expansion, Gram(), gram_frame
  )

  pauli = PauliSpace(:jet_completion_pauli)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  recursive_frame = DissipativeFrame(σx, σy, σz)
  recursive_expansion = floquet_expansion(
    cos(ω * t) * σx, ω, t, VanVleck(), 3; channels=(collapse(σz),)
  )

  JET.@test_opt target_modules=(FloquetExpansions,) positive_completion(
    recursive_expansion, Gram(), recursive_frame
  )

  spectral_frame = DissipativeFrame(σz, σy)
  JET.@test_opt target_modules=(FloquetExpansions,) positive_completion(
    recursive_expansion, Spectral(), spectral_frame
  )

  completion = positive_completion(gram_expansion, Gram(), gram_frame)
  JET.@test_opt target_modules=(FloquetExpansions,) hamiltonian(completion)
  JET.@test_opt target_modules=(FloquetExpansions,) kossakowski_component(completion, 0)
end

@testset "native CP-HFE optimizer stability" begin
  space = NLevelSpace(:jet_cp_hfe, 2)
  σx = Transition(space, :σ, 1, 2) + Transition(space, :σ, 2, 1)
  σz = Transition(space, :σ, 1, 1) - Transition(space, :σ, 2, 2)
  σminus = Transition(space, :σ, 1, 2)
  @variables ω_cp::Real t_cp::Real γ_cp::Real

  H = σz + σx * SQAJet.expim(-ω_cp * t_cp) + σx * SQAJet.expim(ω_cp * t_cp)
  seed = only(FloquetExpansions.physical_amplitude_seeds((jump(σminus, γ_cp),), ω_cp, t_cp))
  coherent = floquet_expansion(harmonics(H, ω_cp, t_cp), VanVleck(), 2)
  kicks = getfield(coherent, :kick_components)
  transported = FloquetExpansions.transport_amplitude_series(seed, kicks, 2)

  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.reattach(kicks[1], 1)

  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.transport_amplitude_series(
    seed, kicks, 2
  )

  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.reconstruct_cp_amplitude_channels([
    transported
  ])

  rows = FloquetExpansions.reconstruct_cp_amplitude_channels([transported])
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.reconstruct_cp_effective_generator(
    coherent, rows
  )
end

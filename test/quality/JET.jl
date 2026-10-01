using Test
using FloquetExpansions
using JET: JET
using Symbolics: @variables

@testset "JET" begin
  JET.test_package(FloquetExpansions; target_modules=(FloquetExpansions,))
end

@testset "expansion algorithm optimizer stability" begin
  pauli = PauliSpace(:jet_expansion_algorithm)
  σx = Pauli(pauli, :sigma, 1)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω_bf::Real
  generator = PeriodicGenerator(Dict(0 => 1 * σz, 1 => 1 * σx, -1 => 1 * σx), ω_bf)

  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    generator, VanVleck(), 3
  )
  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    generator, VanVleck(; algorithm=BlochFeshbach()), 3
  )

  liouvillian_generator = PeriodicGenerator(
    Dict(0 => hamiltonian_action(σz), 1 => dissipator(σx), -1 => hamiltonian_action(σx)),
    ω_bf,
  )
  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    liouvillian_generator, VanVleck(), 2
  )
  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    liouvillian_generator, VanVleck(; algorithm=BlochFeshbach()), 2
  )
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

@testset "native graded channel optimizer stability" begin
  T = ComplexF64
  channels = [
    FloquetExpansions.GradedChannel{T}(0, [T[1, 0, 0], T[0, 1, 0], T[0, 0, 1]]),
    FloquetExpansions.GradedChannel{T}(1, [T[0, 1, 1]]),
  ]
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.active_channels(
    channels, 2, 3
  )
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.known_gram(
    channels, 2, 3
  )
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.gram_coefficient(
    channels, 2, 3
  )
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.store_births!(
    copy(channels), zeros(T, 3, 1), 2
  )
end

@testset "native static slot optimizer stability" begin
  T = ComplexF64
  residual = T[2 1im 0; -1im 1 0; 0 0 0.5]
  active = reshape(T[1, 0, 0], :, 1)
  images = [Matrix{T}(FloquetExpansions.LinearAlgebra.I, 3, 3)]
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.native_static_solve(
    residual, zeros(T, 3, 3), active, images, 1e-8
  )
end

@testset "native static step optimizer stability" begin
  rep = FloquetExpansions.DenseLiouvilleRepresentation(2)
  T = ComplexF64
  σm = T[0 0; 1 0]
  L0 = FloquetExpansions.native_gksl(rep, T[0.5 0; 0 -0.5], zeros(T, 3, 3))
  active = reshape(T[0.5, 0.5im, 0], :, 1)
  residual = FloquetExpansions.native_gksl(rep, zeros(T, 2, 2), active * active')
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.native_static_step(
    rep,
    FloquetExpansions.NoHomologicalInverse(),
    L0,
    residual,
    zeros(T, 3, 3),
    active,
    1e-8,
  )
end

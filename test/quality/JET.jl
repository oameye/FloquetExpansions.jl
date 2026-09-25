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
    liouvillian_generator,
    VanVleck(; algorithm=BlochFeshbach(; complete_positive=Val(false))),
    2,
  )
end

@testset "CP algorithm policy optimizer stability" begin
  JET.@test_opt target_modules=(FloquetExpansions,) HoriDeprit(;
    complete_positive=Val(false)
  )
  JET.@test_opt target_modules=(FloquetExpansions,) BlochFeshbach(;
    complete_positive=Val(false)
  )
  JET.@test_opt target_modules=(FloquetExpansions,) HoriDeprit(;
    complete_positive=Val(true)
  )
  JET.@test_opt target_modules=(FloquetExpansions,) BlochFeshbach(;
    complete_positive=Val(true)
  )

  pauli = PauliSpace(:jet_cp_algorithm_policy)
  σx = Pauli(pauli, :sigma, 1)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω_cp_algorithm_policy::Real
  hamiltonian_generator = PeriodicGenerator(
    Dict(0 => 1 * σz, 1 => 1 * σx, -1 => 1 * σx), ω_cp_algorithm_policy
  )
  liouvillian_generator = PeriodicGenerator(
    Dict(0 => hamiltonian_action(σz) + dissipator(σx)), ω_cp_algorithm_policy
  )

  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    hamiltonian_generator,
    VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(true))),
    2,
  )
  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    hamiltonian_generator,
    VanVleck(; algorithm=BlochFeshbach(; complete_positive=Val(true))),
    2,
  )
  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    liouvillian_generator, VanVleck(; algorithm=HoriDeprit()), 2
  )
  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    liouvillian_generator, VanVleck(; algorithm=BlochFeshbach()), 2
  )
  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    liouvillian_generator,
    VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))),
    2,
  )
  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    liouvillian_generator,
    VanVleck(; algorithm=BlochFeshbach(; complete_positive=Val(false))),
    2,
  )
end

@testset "completion optimizer stability" begin
  fock = FockSpace(:jet_completion_fock)
  a = Destroy(fock, :a)
  @variables ω::Real t::Real
  gram_frame = DissipativeFrame(a, a^2)
  gram_generator = liouvillian(0 * a; channels=(collapse(a + a^2), collapse(a + im * a^2)))
  gram_expansion = floquet_expansion(
    gram_generator,
    ω,
    t,
    VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))),
    1,
  )

  JET.@test_opt target_modules=(FloquetExpansions,) positive_completion(
    gram_expansion, Gram(), gram_frame
  )

  pauli = PauliSpace(:jet_completion_pauli)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  recursive_frame = DissipativeFrame(σx, σy, σz)
  recursive_expansion = floquet_expansion(
    cos(ω * t) * σx,
    ω,
    t,
    VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))),
    3;
    channels=(collapse(σz),),
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

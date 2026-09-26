using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA_EA = SecondQuantizedAlgebra

struct UnsupportedExpansionAlgorithm <: ExpansionAlgorithm end

function ea_vanishes(value)
  return iszero(SQA_EA.simplify(value))
end

function ea_vanishes(generator::PeriodicGenerator)
  return iszero(SQA_EA.simplify(generator))
end

@testset "Van Vleck algorithm selectors" begin
  default = VanVleck()
  explicit_hori_deprit = VanVleck(; algorithm=HoriDeprit())
  bloch_feshbach = VanVleck(; algorithm=BlochFeshbach())

  @test default == explicit_hori_deprit
  @test default != bloch_feshbach
  @test typeof(default) === typeof(explicit_hori_deprit)
  @test typeof(default) !== typeof(bloch_feshbach)
end

@testset "Bloch Feshbach matches Hamiltonian Hori Deprit Van Vleck" begin
  space = PauliSpace(:expansion_algorithm_hamiltonian)
  σx = Pauli(space, :σ, 1)
  σy = Pauli(space, :σ, 2)
  σz = Pauli(space, :σ, 3)
  @variables ω_ea::Real

  H1 = σx + im * σy
  H2 = 2 * σx - im * σz
  H = PeriodicGenerator(Dict(0 => 1 * σz, 1 => H1, -1 => H1', 2 => H2, -2 => H2'), ω_ea)
  order = 4

  default = floquet_expansion(H, VanVleck(), order)
  explicit_hori_deprit = floquet_expansion(H, VanVleck(; algorithm=HoriDeprit()), order)
  bloch_feshbach = floquet_expansion(H, VanVleck(; algorithm=BlochFeshbach()), order)

  for n in 0:(order - 1)
    @test ea_vanishes(
      effective_component(default, n) - effective_component(explicit_hori_deprit, n)
    )
    @test ea_vanishes(
      effective_component(bloch_feshbach, n) - effective_component(default, n)
    )
  end
  for n in 1:(order - 1)
    @test ea_vanishes(micromotion(default, n) - micromotion(explicit_hori_deprit, n))
    @test ea_vanishes(micromotion(bloch_feshbach, n) - micromotion(default, n))
  end

  leading = floquet_expansion(H, VanVleck(; algorithm=BlochFeshbach()), 1)
  @test ea_vanishes(effective_component(leading, 0) - H[0])
  @test iszero(micromotion(leading))

  unsupported = VanVleck(; algorithm=UnsupportedExpansionAlgorithm())
  @test_throws ArgumentError floquet_expansion(H, unsupported, 2)
  @test_throws ArgumentError floquet_expansion(H, VanVleck(; algorithm=BlochFeshbach()), 0)
end

@testset "Bloch Feshbach Hamiltonian leading correction" begin
  space = PauliSpace(:expansion_algorithm_hamiltonian_reference)
  σx = Pauli(space, :σ, 1)
  σy = Pauli(space, :σ, 2)
  σz = Pauli(space, :σ, 3)
  @variables ω_ea_ref::Real

  H1 = σx + im * σy
  H2 = 2 * σx + σz
  H = PeriodicGenerator(Dict(0 => 1 * σz, 1 => H1, -1 => H1', 2 => H2, -2 => H2), ω_ea_ref)
  nonzero_harmonics = filter(!iszero, collect(keys(H)))

  expected_effective = zero(H[0])
  expected_kick_components = Dict{Int,typeof(H[0])}()
  for harmonic in nonzero_harmonics
    expected_effective += (1 // harmonic) * H[-harmonic] * H[harmonic]
    expected_kick_components[harmonic] = (im // harmonic) * H[harmonic]
  end
  expected_kick = PeriodicGenerator(expected_kick_components, ω_ea_ref)

  expansion = floquet_expansion(H, VanVleck(; algorithm=BlochFeshbach()), 2)
  @test ea_vanishes(effective_component(expansion, 1) - ω_ea_ref^(-1) * expected_effective)
  @test ea_vanishes(micromotion(expansion, 1) - ω_ea_ref^(-1) * expected_kick)
end

@testset "Bloch Feshbach matches Liouvillian Hori Deprit" begin
  space = PauliSpace(:expansion_algorithm_liouvillian)
  σx = Pauli(space, :σ, 1)
  σy = Pauli(space, :σ, 2)
  σz = Pauli(space, :σ, 3)
  @variables ω_ea_map::Real

  L = PeriodicGenerator(
    Dict(
      0 => hamiltonian_action(σz) + dissipator(σx),
      1 => hamiltonian_action(σx + σy),
      -1 => dissipator(σy + σz),
    ),
    ω_ea_map,
  )
  order = 3

  hori_deprit = floquet_expansion(
    L, VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))), order
  )
  bloch_feshbach = floquet_expansion(
    L, VanVleck(; algorithm=BlochFeshbach(; complete_positive=Val(false))), order
  )

  for n in 0:(order - 1)
    @test ea_vanishes(
      effective_component(bloch_feshbach, n) - effective_component(hori_deprit, n)
    )
  end
  for n in 1:(order - 1)
    @test ea_vanishes(micromotion(bloch_feshbach, n) - micromotion(hori_deprit, n))
  end
end

@testset "Bloch Feshbach Liouvillian leading correction" begin
  space = PauliSpace(:expansion_algorithm_liouvillian_reference)
  σx = Pauli(space, :σ, 1)
  σy = Pauli(space, :σ, 2)
  σz = Pauli(space, :σ, 3)
  @variables ω_ea_map_ref::Real

  L = PeriodicGenerator(
    Dict(
      0 => hamiltonian_action(σz) + dissipator(σx),
      1 => hamiltonian_action(σx),
      -1 => dissipator(σy),
      2 => hamiltonian_action(σy),
      -2 => dissipator(σz),
    ),
    ω_ea_map_ref,
  )
  nonzero_harmonics = filter(!iszero, collect(keys(L)))

  expected_effective = zero(L[0])
  expected_kick_components = Dict{Int,typeof(L[0])}()
  for harmonic in nonzero_harmonics
    expected_effective += (im // harmonic) * compose(L[-harmonic], L[harmonic])
    expected_kick_components[harmonic] = (im // harmonic) * L[harmonic]
  end
  expected_kick = PeriodicGenerator(expected_kick_components, ω_ea_map_ref)

  expansion = floquet_expansion(
    L, VanVleck(; algorithm=BlochFeshbach(; complete_positive=Val(false))), 2
  )
  @test ea_vanishes(
    effective_component(expansion, 1) - ω_ea_map_ref^(-1) * expected_effective
  )
  @test ea_vanishes(micromotion(expansion, 1) - ω_ea_map_ref^(-1) * expected_kick)
end

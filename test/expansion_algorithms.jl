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
  raw_hori_deprit = VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false)))
  raw_bloch_feshbach = VanVleck(; algorithm=BlochFeshbach(; complete_positive=Val(false)))

  @test default == explicit_hori_deprit
  @test default.algorithm.algorithm isa HoriDeprit
  @test explicit_hori_deprit.algorithm.algorithm isa HoriDeprit
  @test bloch_feshbach.algorithm.algorithm isa BlochFeshbach
  @test raw_hori_deprit.algorithm isa HoriDeprit
  @test raw_bloch_feshbach.algorithm isa BlochFeshbach
  @test typeof(default) === typeof(explicit_hori_deprit)
  @test typeof(default) !== typeof(bloch_feshbach)
  @test typeof(default) !== typeof(raw_hori_deprit)
  @test typeof(bloch_feshbach) !== typeof(raw_bloch_feshbach)
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

  @test getfield(default, :effective_components) ==
    getfield(explicit_hori_deprit, :effective_components)
  @test getfield(default, :kick_components) ==
    getfield(explicit_hori_deprit, :kick_components)

  for n in 0:(order - 1)
    @test ea_vanishes(
      effective_component(bloch_feshbach, n) - effective_component(default, n)
    )
  end
  for n in 1:(order - 1)
    @test ea_vanishes(micromotion(bloch_feshbach, n) - micromotion(default, n))
  end

  leading = floquet_expansion(H, VanVleck(; algorithm=BlochFeshbach()), 1)
  @test ea_vanishes(effective_component(leading, 0) - H[0])
  @test iszero(micromotion(leading))

  unsupported = VanVleck(; algorithm=UnsupportedExpansionAlgorithm())
  @test_throws ArgumentError floquet_expansion(H, unsupported, 2)
end

@testset "Bloch Feshbach matches generic Liouvillian Hori Deprit Van Vleck" begin
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

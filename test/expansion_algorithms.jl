using Test
using FloquetExpansions
using LinearAlgebra: norm
using Random: MersenneTwister
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables

include(joinpath(@__DIR__, "helpers", "shared.jl"))

const SQA_EA = SecondQuantizedAlgebra

struct UnsupportedExpansionAlgorithm <: ExpansionAlgorithm end

function ea_vanishes(value)
  return iszero(SQA_EA.simplify(value))
end

function ea_vanishes(generator::PeriodicGenerator)
  return iszero(SQA_EA.simplify(generator))
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
  order = 6

  hori_deprit = floquet_expansion(H, VanVleck(), order)
  bloch_feshbach = floquet_expansion(H, VanVleck(; algorithm=BlochFeshbach()), order)

  for n in 0:(order - 1)
    @test ea_vanishes(
      effective_component(bloch_feshbach, n) - effective_component(hori_deprit, n)
    )
  end
  for n in 1:(order - 1)
    @test ea_vanishes(micromotion(bloch_feshbach, n) - micromotion(hori_deprit, n))
    @test ea_vanishes(time_average(micromotion(bloch_feshbach, n)))
  end

  leading = floquet_expansion(H, VanVleck(; algorithm=BlochFeshbach()), 1)
  @test ea_vanishes(effective_component(leading, 0) - H[0])
  @test iszero(micromotion(leading))

  unsupported = VanVleck(; algorithm=UnsupportedExpansionAlgorithm())
  @test_throws ArgumentError floquet_expansion(H, unsupported, 2)
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
  expected_micromotion_components = Dict{Int,typeof(H[0])}()
  for harmonic in nonzero_harmonics
    expected_effective += (1 // harmonic) * H[-harmonic] * H[harmonic]
    expected_micromotion_components[harmonic] = (im // harmonic) * H[harmonic]
  end
  expected_micromotion = PeriodicGenerator(expected_micromotion_components, ω_ea_ref)

  expansion = floquet_expansion(H, VanVleck(; algorithm=BlochFeshbach()), 2)
  @test ea_vanishes(effective_component(expansion, 1) - ω_ea_ref^(-1) * expected_effective)
  @test ea_vanishes(micromotion(expansion, 1) - ω_ea_ref^(-1) * expected_micromotion)
end

@testset "Van Vleck second order matches des Cloizeaux" begin
  space = FockSpace(:expansion_algorithm_des_cloizeaux)
  a = Destroy(space, :a)
  @variables ω_ea_dc::Real

  H1 = a + a' * a'
  H = PeriodicGenerator(Dict(0 => a' * a, 1 => H1, -1 => H1'), ω_ea_dc)
  nonzero_harmonics = filter(!iszero, collect(keys(H)))

  brillouin_wigner = zero(H[0])
  metric = zero(H[0])
  for m in nonzero_harmonics
    for n in nonzero_harmonics
      brillouin_wigner += (1 // (m * n)) * H[-m] * H[m - n] * H[n]
    end
    brillouin_wigner -= (1 // m^2) * H[-m] * H[m] * H[0]
    metric += (1 // m^2) * H[-m] * H[m]
  end
  metric_correction = (1 // 2) * (metric * H[0] - H[0] * metric)
  des_cloizeaux = brillouin_wigner + metric_correction
  @test !ea_vanishes(metric_correction)

  for algorithm in (HoriDeprit(), BlochFeshbach())
    expansion = floquet_expansion(H, VanVleck(; algorithm), 3)
    @test ea_vanishes(effective_component(expansion, 2) - ω_ea_dc^(-2) * des_cloizeaux)
  end
end

@testset "Bloch Feshbach matches Liouvillian Hori Deprit" begin
  space = PauliSpace(:expansion_algorithm_liouvillian)
  σx = Pauli(space, :σ, 1)
  σy = Pauli(space, :σ, 2)
  σz = Pauli(space, :σ, 3)
  @variables ω_ea_map::Real

  symmetric = PeriodicGenerator(
    Dict(
      0 => hamiltonian_action(σz) + dissipator(σx),
      1 => hamiltonian_action(σx + σy),
      -1 => dissipator(σy + σz),
    ),
    ω_ea_map,
  )
  resonant_triple = PeriodicGenerator(
    Dict(
      0 => hamiltonian_action(σz) + dissipator(σx),
      1 => hamiltonian_action(σx + σy),
      2 => dissipator(σy),
      -3 => hamiltonian_action(σz) + dissipator(σx + σz),
    ),
    ω_ea_map,
  )
  order = 5

  for L in (symmetric, resonant_triple)
    hori_deprit = floquet_expansion(L, VanVleck(; algorithm=HoriDeprit()), order)
    bloch_feshbach = floquet_expansion(L, VanVleck(; algorithm=BlochFeshbach()), order)

    for n in 0:(order - 1)
      @test ea_vanishes(
        effective_component(bloch_feshbach, n) - effective_component(hori_deprit, n)
      )
    end
    for n in 1:(order - 1)
      @test ea_vanishes(micromotion(bloch_feshbach, n) - micromotion(hori_deprit, n))
      @test ea_vanishes(time_average(micromotion(bloch_feshbach, n)))
    end
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
  expected_micromotion_components = Dict{Int,typeof(L[0])}()
  for harmonic in nonzero_harmonics
    expected_effective += (im // harmonic) * compose(L[-harmonic], L[harmonic])
    expected_micromotion_components[harmonic] = (im // harmonic) * L[harmonic]
  end
  expected_micromotion = PeriodicGenerator(expected_micromotion_components, ω_ea_map_ref)

  expansion = floquet_expansion(L, VanVleck(; algorithm=BlochFeshbach()), 2)
  @test ea_vanishes(
    effective_component(expansion, 1) - ω_ea_map_ref^(-1) * expected_effective
  )
  @test ea_vanishes(micromotion(expansion, 1) - ω_ea_map_ref^(-1) * expected_micromotion)
end

@testset "Bloch Feshbach expansion completes like Hori Deprit" begin
  space = PauliSpace(:expansion_algorithm_completion)
  σx = Pauli(space, :σ, 1)
  σy = Pauli(space, :σ, 2)
  σz = Pauli(space, :σ, 3)
  σminus = (1 // 2) * (σx - im * σy)
  @variables ω_ea_cp::Real t_ea_cp::Real E_ea_cp::Real γ_ea_cp::Real

  H = (1 // 2) * σz + E_ea_cp * cos(ω_ea_cp * t_ea_cp) * σx
  channels = (jump(σminus, γ_ea_cp),)
  frame = DissipativeFrame(σx, σy, σz)

  hori_deprit, bloch_feshbach = map((HoriDeprit(), BlochFeshbach())) do algorithm
    raw = floquet_expansion(
      H, ω_ea_cp, t_ea_cp, VanVleck(; algorithm), 2; channels=channels
    )
    positive_completion(raw, Gram(), frame)
  end

  @test ea_vanishes(effective_generator(bloch_feshbach) - effective_generator(hori_deprit))
end

function ea_numeric(coefficient, substitutions)
  value = Symbolics.substitute(SQA_EA.to_num(coefficient), substitutions)
  return complex(
    Float64(Symbolics.value(real(value))), Float64(Symbolics.value(imag(value)))
  )
end

function ea_superoperator(L::Liouvillian, d::Int, substitutions)
  S = zeros(ComplexF64, d^2, d^2)
  for (left, right, coefficient) in terms(L)
    left_matrix = tomatrix(left, d, substitutions)
    right_matrix = tomatrix(right, d, substitutions)
    S .+=
      ea_numeric(coefficient, substitutions) .* kron(transpose(right_matrix), left_matrix)
  end
  return S
end

@testset "Coherent Liouvillian expansion matches the Hamiltonian expansion" begin
  d = 3
  atom = NLevelSpace(:expansion_algorithm_coherent, d)
  @variables ω_ea_coh::Real
  H = random_drive(MersenneTwister(0x0AC1E), atom, d, 2, ω_ea_coh)
  L = PeriodicGenerator(Dict(m => hamiltonian_action(H[m]) for m in keys(H)), ω_ea_coh)
  substitutions = Dict(ω_ea_coh => 7.0)
  order = 4

  superoperator(X) = ea_superoperator(X, d, substitutions)
  function coherent_superoperator(X)
    return ea_superoperator(hamiltonian_action(X), d, substitutions)
  end

  for algorithm in (HoriDeprit(), BlochFeshbach())
    hamiltonian = floquet_expansion(H, VanVleck(; algorithm), order)
    liouvillian = floquet_expansion(L, VanVleck(; algorithm), order)

    for n in 0:(order - 1)
      expected = coherent_superoperator(effective_component(hamiltonian, n))
      actual = superoperator(effective_component(liouvillian, n))
      @test norm(actual - expected) <= 1e-10 * max(1, norm(expected))
    end
    for n in 1:(order - 1)
      coherent = micromotion(hamiltonian, n)
      dissipative = micromotion(liouvillian, n)
      for m in union(keys(coherent), keys(dissipative))
        expected = coherent_superoperator(coherent[m])
        actual = superoperator(dissipative[m])
        @test norm(actual - expected) <= 1e-10 * max(1, norm(expected))
      end
    end
  end
end

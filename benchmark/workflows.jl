using BenchmarkTools: @benchmarkable
using FloquetExpansions
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

function kerr_parametric_oscillator()
  space = FockSpace(:cavity)
  a = Destroy(space, :a)
  @variables ω::Real t::Real δ::Real K::Real A::Real B::Real

  H =
    δ * a' * a +
    (K / 2) * a' * a' * a * a +
    A * cos(ω * t) * (a + a') +
    B * cos(2 * ω * t) * (a + a')
  return H, ω, t
end

function driven_qubit()
  space = NLevelSpace(:qubit, 2)
  σ11 = Transition(space, :σ, 1, 1)
  σ22 = Transition(space, :σ, 2, 2)
  σ12 = Transition(space, :σ, 1, 2)
  σ21 = Transition(space, :σ, 2, 1)
  σz = σ11 - σ22
  σx = σ12 + σ21
  @variables ω::Real t::Real Δ::Real A::Real

  H = (Δ / 2) * σz + A * cos(ω * t) * σx
  return H, ω, t
end

function driven_dissipative_qubit()
  pauli = PauliSpace(:dissipative_qubit_benchmark)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω::Real t::Real Δ::Real A::Real γ::Real

  H = (Δ / 2) * σz + A * cos(ω * t) * σx
  L = hamiltonian_action(H) + γ * dissipator((σx - im * σy) / 2)
  return L, ω, t
end

function gram_completion_workload()
  space = FockSpace(:completion_benchmark)
  a = Destroy(space, :a)
  @variables ω::Real t::Real
  frame = JumpBasis(a, a^2)
  generator = SQA.SuperOperator(dissipator(a + a^2) + dissipator(a + im * a^2))
  expansion = floquet_expansion(generator, ω, t, VanVleck(), 1)
  return expansion, frame
end

function recursive_gram_workload()
  pauli = PauliSpace(:recursive_completion_benchmark)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω::Real t::Real Ω::Real
  frame = JumpBasis(σx, σy, σz)
  expansion = floquet_expansion(
    hamiltonian_action(Ω * cos(ω * t) * σx) + dissipator(σz), ω, t, VanVleck(), 3
  )
  return expansion, frame
end

function spectral_completion_workload()
  pauli = PauliSpace(:spectral_completion_benchmark)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω::Real t::Real Ω::Real
  frame = JumpBasis(σz, σy)
  expansion = floquet_expansion(
    hamiltonian_action(Ω * cos(ω * t) * σx) + dissipator(σz), ω, t, VanVleck(), 3
  )
  return expansion, frame
end

function benchmark_fourier_expansion!(suite)
  kpo, kpo_ω, kpo_t = kerr_parametric_oscillator()
  qubit, qubit_ω, qubit_t = driven_qubit()

  suite["Fourier Expansion"]["Kerr parametric oscillator"]["symbolic input"] = @benchmarkable harmonics(
    $kpo, $kpo_ω, $kpo_t
  )
  suite["Fourier Expansion"]["Driven qubit"]["symbolic input"] = @benchmarkable harmonics(
    $qubit, $qubit_ω, $qubit_t
  )
  return nothing
end

function benchmark_floquet_expansion!(suite)
  kpo, kpo_ω, kpo_t = kerr_parametric_oscillator()
  qubit, qubit_ω, qubit_t = driven_qubit()

  for order in 1:3
    suite["Floquet Expansion"]["Kerr parametric oscillator"]["order $order"] = @benchmarkable begin
      vv = floquet_expansion($kpo, $kpo_ω, $kpo_t, VanVleck(), $order)
      effective_generator(vv), micromotion(vv)
    end
    suite["Floquet Expansion"]["Driven qubit"]["order $order"] = @benchmarkable begin
      vv = floquet_expansion($qubit, $qubit_ω, $qubit_t, VanVleck(), $order)
      effective_generator(vv), micromotion(vv)
    end
  end
  return nothing
end

function benchmark_expansion_algorithms!(suite)
  kpo, kpo_ω, kpo_t = kerr_parametric_oscillator()
  qubit, qubit_ω, qubit_t = driven_qubit()
  open_qubit, open_ω, open_t = driven_dissipative_qubit()
  algorithms = (("Hori-Deprit", HoriDeprit()), ("Bloch-Feshbach", BlochFeshbach()))

  for (name, algorithm) in algorithms, order in (3, 5)
    gauge = VanVleck(; algorithm)
    suite["Expansion Algorithm"]["Kerr parametric oscillator"]["$name order $order"] = @benchmarkable floquet_expansion(
      $kpo, $kpo_ω, $kpo_t, $gauge, $order
    )
    suite["Expansion Algorithm"]["Driven qubit"]["$name order $order"] = @benchmarkable floquet_expansion(
      $qubit, $qubit_ω, $qubit_t, $gauge, $order
    )
    suite["Expansion Algorithm"]["Driven dissipative qubit"]["$name order $order"] = @benchmarkable floquet_expansion(
      $open_qubit, $open_ω, $open_t, $gauge, $order
    )
  end
  return nothing
end

function benchmark_positive_completion!(suite)
  expansion, frame = gram_completion_workload()
  recursive_expansion, recursive_frame = recursive_gram_workload()
  spectral_expansion, spectral_frame = spectral_completion_workload()

  suite["Positive Completion"]["Full-rank bosonic"]["Gram fixed frame"] = @benchmarkable positive_completion(
    $expansion, Gram(), $frame
  )
  suite["Positive Completion"]["Full-rank bosonic"]["Gram automatic frame"] = @benchmarkable positive_completion(
    $expansion, Gram()
  )
  suite["Positive Completion"]["Recursive dark onset"]["Gram fixed frame"] = @benchmarkable positive_completion(
    $recursive_expansion, Gram(), $recursive_frame
  )
  suite["Positive Completion"]["Driven qubit"]["Gram fixed frame"] = @benchmarkable positive_completion(
    $spectral_expansion, Gram(), $spectral_frame
  )
  suite["Positive Completion"]["Driven qubit"]["Spectral fixed frame"] = @benchmarkable positive_completion(
    $spectral_expansion, Spectral(), $spectral_frame
  )
  return nothing
end

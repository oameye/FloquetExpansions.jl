using BenchmarkTools: @benchmarkable
using FloquetExpansions
using Symbolics: @variables

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

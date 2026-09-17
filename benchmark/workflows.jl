using BenchmarkTools: @benchmarkable
using FloquetExpansions
using Symbolics: @variables

const FE = FloquetExpansions

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

function gram_completion_workload()
  space = FockSpace(:completion_benchmark)
  a = Destroy(space, :a)
  @variables ω::Real t::Real
  frame = DissipativeFrame(a, a^2)
  generator = liouvillian(0 * a; channels=(collapse(a + a^2), collapse(a + im * a^2)))
  expansion = floquet_expansion(generator, ω, t, VanVleck(), 1)
  return expansion, frame
end

function recursive_gram_workload()
  pauli = PauliSpace(:recursive_completion_benchmark)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω::Real t::Real Ω::Real
  frame = DissipativeFrame(σx, σy, σz)
  expansion = floquet_expansion(
    Ω * cos(ω * t) * σx, ω, t, VanVleck(), 3; channels=(collapse(σz),)
  )
  return expansion, frame
end

function spectral_completion_workload()
  pauli = PauliSpace(:spectral_completion_benchmark)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω::Real t::Real Ω::Real
  frame = DissipativeFrame(σz, σy)
  expansion = floquet_expansion(
    Ω * cos(ω * t) * σx, ω, t, VanVleck(), 3; channels=(collapse(σz),)
  )
  return expansion, frame
end

function cp_hfe_qubit_workload()
  pauli = PauliSpace(:cp_hfe_benchmark_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  bright = σy + σz
  dark = σy - σz
  @variables ω_cp_bench::Real t_cp_bench::Real Ω_cp_bench::Real

  ω = ω_cp_bench
  t = t_cp_bench
  H = Ω_cp_bench * cos(ω * t) * σx
  channels = (collapse(bright),)
  frame = DissipativeFrame(bright, dark)
  return H, ω, t, channels, frame
end

function cp_hfe_bosonic_workload()
  fock = FockSpace(:cp_hfe_benchmark_kerr)
  a = Destroy(fock, :a)
  @variables ω_cp_boson::Real t_cp_boson::Real Δ_cp_boson::Real K_cp_boson::Real ε_cp_boson::Real

  ω = ω_cp_boson
  t = t_cp_boson
  number = a' * a
  H =
    Δ_cp_boson * number +
    K_cp_boson * a'^2 * a^2 +
    ε_cp_boson * cos(ω * t) * (a + a')
  channels = (collapse(a^2),)
  return H, ω, t, channels
end

function cp_hfe_stage_context(H, ω, t, order, channels)
  coherent = floquet_expansion(H, ω, t, VanVleck(), order)
  kicks = getfield(coherent, :kick_components)
  seeds = FE.physical_amplitude_seeds(channels, ω, t)
  amplitudes = [FE.transport_amplitude_series(seed, kicks, order) for seed in seeds]
  reconstructed_channels = FE.reconstruct_cp_amplitude_channels(amplitudes)
  return coherent, kicks, seeds, amplitudes, reconstructed_channels
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

function benchmark_positive_completion!(suite)
  expansion, frame = gram_completion_workload()
  recursive_expansion, recursive_frame = recursive_gram_workload()
  spectral_expansion, spectral_frame = spectral_completion_workload()

  suite["Positive Completion"]["Full-rank bosonic"]["fixed frame"]["Gram"] = @benchmarkable positive_completion(
    $expansion, Gram(), $frame
  )
  suite["Positive Completion"]["Full-rank bosonic"]["automatic frame"]["Gram"] = @benchmarkable positive_completion(
    $expansion, Gram()
  )
  suite["Positive Completion"]["Recursive dark onset"]["fixed frame"]["Gram"] = @benchmarkable positive_completion(
    $recursive_expansion, Gram(), $recursive_frame
  )
  suite["Positive Completion"]["Driven qubit"]["fixed frame"]["Gram"] = @benchmarkable positive_completion(
    $spectral_expansion, Gram(), $spectral_frame
  )
  suite["Positive Completion"]["Driven qubit"]["fixed frame"]["Spectral"] = @benchmarkable positive_completion(
    $spectral_expansion, Spectral(), $spectral_frame
  )
  return nothing
end

function benchmark_cp_hfe!(suite)
  qubit_H, qubit_ω, qubit_t, qubit_channels, qubit_frame = cp_hfe_qubit_workload()
  boson_H, boson_ω, boson_t, boson_channels = cp_hfe_bosonic_workload()

  for order in 1:3
    suite["CP HFE"]["Driven qubit"]["end to end"]["order $order"] = @benchmarkable FE.cp_hfe_reconstruction(
      $qubit_H, $qubit_ω, $qubit_t, $order, $qubit_channels
    )
    suite["CP HFE"]["Driven Kerr resonator"]["end to end"]["order $order"] = @benchmarkable FE.cp_hfe_reconstruction(
      $boson_H, $boson_ω, $boson_t, $order, $boson_channels
    )
  end

  order = 3
  qubit_coherent, qubit_kicks, qubit_seeds, qubit_amplitudes, qubit_cp_channels =
    cp_hfe_stage_context(qubit_H, qubit_ω, qubit_t, order, qubit_channels)
  boson_coherent, boson_kicks, boson_seeds, boson_amplitudes, boson_cp_channels =
    cp_hfe_stage_context(boson_H, boson_ω, boson_t, order, boson_channels)

  suite["CP HFE"]["Driven qubit"]["stages"]["coherent solve"] = @benchmarkable floquet_expansion(
    $qubit_H, $qubit_ω, $qubit_t, VanVleck(), $order
  )
  suite["CP HFE"]["Driven qubit"]["stages"]["physical amplitude seeds"] = @benchmarkable FE.physical_amplitude_seeds(
    $qubit_channels, $qubit_ω, $qubit_t
  )
  suite["CP HFE"]["Driven qubit"]["stages"]["amplitude transport"] = @benchmarkable [
    FE.transport_amplitude_series(seed, $qubit_kicks, $order) for seed in $qubit_seeds
  ]
  suite["CP HFE"]["Driven qubit"]["stages"]["Kraus rows"] = @benchmarkable FE.reconstruct_cp_amplitude_channels(
    $qubit_amplitudes
  )
  suite["CP HFE"]["Driven qubit"]["stages"]["GKSL reconstruction"] = @benchmarkable FE.reconstruct_cp_effective_generator(
    $qubit_coherent, $qubit_cp_channels
  )

  suite["CP HFE"]["Driven Kerr resonator"]["stages"]["coherent solve"] = @benchmarkable floquet_expansion(
    $boson_H, $boson_ω, $boson_t, VanVleck(), $order
  )
  suite["CP HFE"]["Driven Kerr resonator"]["stages"]["physical amplitude seeds"] = @benchmarkable FE.physical_amplitude_seeds(
    $boson_channels, $boson_ω, $boson_t
  )
  suite["CP HFE"]["Driven Kerr resonator"]["stages"]["amplitude transport"] = @benchmarkable [
    FE.transport_amplitude_series(seed, $boson_kicks, $order) for seed in $boson_seeds
  ]
  suite["CP HFE"]["Driven Kerr resonator"]["stages"]["Kraus rows"] = @benchmarkable FE.reconstruct_cp_amplitude_channels(
    $boson_amplitudes
  )
  suite["CP HFE"]["Driven Kerr resonator"]["stages"]["GKSL reconstruction"] = @benchmarkable FE.reconstruct_cp_effective_generator(
    $boson_coherent, $boson_cp_channels
  )

  suite["CP HFE Comparison"]["Driven qubit"]["order 3"]["native CP HFE"] = @benchmarkable FE.cp_hfe_reconstruction(
    $qubit_H, $qubit_ω, $qubit_t, $order, $qubit_channels
  )
  suite["CP HFE Comparison"]["Driven qubit"]["order 3"]["raw Liouvillian Van Vleck"] = @benchmarkable begin
    raw = floquet_expansion(
      $qubit_H, $qubit_ω, $qubit_t, VanVleck(), $order; channels=$qubit_channels
    )
    effective_generator(raw), micromotion(raw)
  end
  suite["CP HFE Comparison"]["Driven qubit"]["order 3"]["raw VV + Gram"] = @benchmarkable begin
    raw = floquet_expansion(
      $qubit_H, $qubit_ω, $qubit_t, VanVleck(), $order; channels=$qubit_channels
    )
    positive_completion(raw, Gram(), $qubit_frame)
  end
  suite["CP HFE Comparison"]["Driven qubit"]["order 3"]["raw VV + Spectral"] = @benchmarkable begin
    raw = floquet_expansion(
      $qubit_H, $qubit_ω, $qubit_t, VanVleck(), $order; channels=$qubit_channels
    )
    positive_completion(raw, Spectral(), $qubit_frame)
  end
  return nothing
end
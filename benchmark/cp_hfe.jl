const FE = FloquetExpansions

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
  @variables ω_cp_boson::Real t_cp_boson::Real Δ_cp_boson::Real
  @variables K_cp_boson::Real ε_cp_boson::Real

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

function cp_hfe_transport_all(seeds, kicks, order)
  return [FE.transport_amplitude_series(seed, kicks, order) for seed in seeds]
end

function cp_hfe_stage_context(H, ω, t, order, channels)
  coherent = floquet_expansion(H, ω, t, VanVleck(), order)
  kicks = getfield(coherent, :kick_components)
  seeds = FE.physical_amplitude_seeds(channels, ω, t)
  amplitudes = cp_hfe_transport_all(seeds, kicks, order)
  reconstructed = FE.reconstruct_cp_amplitude_channels(amplitudes)
  return coherent, kicks, seeds, amplitudes, reconstructed
end

function benchmark_cp_hfe!(suite)
  qubit_H, qubit_ω, qubit_t, qubit_channels, qubit_frame = cp_hfe_qubit_workload()
  boson_H, boson_ω, boson_t, boson_channels = cp_hfe_bosonic_workload()

  for order in 1:3
    suite["CP HFE"]["Driven qubit"]["end to end"]["order $order"] =
      @benchmarkable FE.cp_hfe_reconstruction(
        $qubit_H, $qubit_ω, $qubit_t, $order, $qubit_channels
      )
    suite["CP HFE"]["Driven Kerr resonator"]["end to end"]["order $order"] =
      @benchmarkable FE.cp_hfe_reconstruction(
        $boson_H, $boson_ω, $boson_t, $order, $boson_channels
      )
  end

  order = 3
  qubit_coherent, qubit_kicks, qubit_seeds, qubit_amplitudes, qubit_rows =
    cp_hfe_stage_context(qubit_H, qubit_ω, qubit_t, order, qubit_channels)
  boson_coherent, boson_kicks, boson_seeds, boson_amplitudes, boson_rows =
    cp_hfe_stage_context(boson_H, boson_ω, boson_t, order, boson_channels)

  suite["CP HFE"]["Driven qubit"]["stages"]["coherent solve"] =
    @benchmarkable floquet_expansion($qubit_H, $qubit_ω, $qubit_t, VanVleck(), $order)
  suite["CP HFE"]["Driven qubit"]["stages"]["physical amplitude seeds"] =
    @benchmarkable FE.physical_amplitude_seeds($qubit_channels, $qubit_ω, $qubit_t)
  suite["CP HFE"]["Driven qubit"]["stages"]["amplitude transport"] =
    @benchmarkable cp_hfe_transport_all($qubit_seeds, $qubit_kicks, $order)
  suite["CP HFE"]["Driven qubit"]["stages"]["Kraus rows"] =
    @benchmarkable FE.reconstruct_cp_amplitude_channels($qubit_amplitudes)
  suite["CP HFE"]["Driven qubit"]["stages"]["GKSL reconstruction"] =
    @benchmarkable FE.reconstruct_cp_effective_generator($qubit_coherent, $qubit_rows)

  suite["CP HFE"]["Driven Kerr resonator"]["stages"]["coherent solve"] =
    @benchmarkable floquet_expansion($boson_H, $boson_ω, $boson_t, VanVleck(), $order)
  suite["CP HFE"]["Driven Kerr resonator"]["stages"]["physical amplitude seeds"] =
    @benchmarkable FE.physical_amplitude_seeds($boson_channels, $boson_ω, $boson_t)
  suite["CP HFE"]["Driven Kerr resonator"]["stages"]["amplitude transport"] =
    @benchmarkable cp_hfe_transport_all($boson_seeds, $boson_kicks, $order)
  suite["CP HFE"]["Driven Kerr resonator"]["stages"]["Kraus rows"] =
    @benchmarkable FE.reconstruct_cp_amplitude_channels($boson_amplitudes)
  suite["CP HFE"]["Driven Kerr resonator"]["stages"]["GKSL reconstruction"] =
    @benchmarkable FE.reconstruct_cp_effective_generator($boson_coherent, $boson_rows)

  suite["CP HFE Comparison"]["Driven qubit"]["order 3"]["native CP HFE"] =
    @benchmarkable FE.cp_hfe_reconstruction(
      $qubit_H, $qubit_ω, $qubit_t, $order, $qubit_channels
    )
  suite["CP HFE Comparison"]["Driven qubit"]["order 3"]["raw Liouvillian Van Vleck"] =
    @benchmarkable begin
      raw = floquet_expansion(
        $qubit_H, $qubit_ω, $qubit_t, VanVleck(), $order; channels=$qubit_channels
      )
      effective_generator(raw), micromotion(raw)
    end
  suite["CP HFE Comparison"]["Driven qubit"]["order 3"]["raw VV + Gram"] =
    @benchmarkable begin
      raw = floquet_expansion(
        $qubit_H, $qubit_ω, $qubit_t, VanVleck(), $order; channels=$qubit_channels
      )
      positive_completion(raw, Gram(), $qubit_frame)
    end
  suite["CP HFE Comparison"]["Driven qubit"]["order 3"]["raw VV + Spectral"] =
    @benchmarkable begin
      raw = floquet_expansion(
        $qubit_H, $qubit_ω, $qubit_t, VanVleck(), $order; channels=$qubit_channels
      )
      positive_completion(raw, Spectral(), $qubit_frame)
    end
  return nothing
end

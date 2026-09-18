mutable struct CoherentWorkProfile
  dressed_generator_nodes::Int
  dressed_kick_nodes::Int
  periodic_commutators::Int
  harmonic_commutator_pairs::Int
  assembly_terms::Int
  assembly_harmonic_inputs::Int
  periodic_simplify_calls::Int
  component_simplify_calls::Int
  simplify_harmonic_inputs::Int
  materialized_harmonics::Int
  peak_harmonics::Int
end

CoherentWorkProfile() = CoherentWorkProfile(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

function observe_harmonics!(profile::CoherentWorkProfile, generator)
  harmonics = length(generator)
  profile.materialized_harmonics += harmonics
  profile.peak_harmonics = max(profile.peak_harmonics, harmonics)
  return generator
end

function record_commutator!(profile::CoherentWorkProfile, left, right)
  profile.periodic_commutators += 1
  profile.harmonic_commutator_pairs += length(left) * length(right)
  return nothing
end

function record_generator_node!(profile, K, dressed_generator, n, j)
  profile.dressed_generator_nodes += 1
  for k in 1:(n - j + 1)
    previous = dressed_generator[FE.triindex(n - k, j - 1)]
    record_commutator!(profile, K[k], previous)
  end
  return nothing
end

function record_kick_node!(profile, K, Kdot, dressed_kick_derivative, n, j)
  profile.dressed_kick_nodes += 1
  for k in 1:(n - j + 1)
    previous = if j == 1
      Kdot[n - k + 1]
    else
      dressed_kick_derivative[FE.triindex(n - k, j - 1)]
    end
    record_commutator!(profile, K[k], previous)
  end
  return nothing
end

function record_assembly!(profile, dressed_generator, dressed_kick_derivative, n)
  for j in 0:n
    profile.assembly_terms += 1
    profile.assembly_harmonic_inputs += length(dressed_generator[FE.triindex(n, j)])
  end
  for j in 1:n
    profile.assembly_terms += 1
    profile.assembly_harmonic_inputs += length(dressed_kick_derivative[FE.triindex(n, j)])
  end
  return nothing
end

function coherent_work_profile(generator, gauge, order)
  nodes = (order * (order + 1)) ÷ 2
  dressed_generator = [zero(generator) for _ in 1:nodes]
  dressed_kick_derivative = [zero(generator) for _ in 1:nodes]
  generator_type = typeof(generator)
  K = generator_type[]
  Kdot = generator_type[]
  effective = typeof(FE.time_average(generator))[]
  profile = CoherentWorkProfile()
  observe_harmonics!(profile, generator)

  for n in 0:(order - 1)
    dressed_generator[FE.triindex(n, 0)] = n == 0 ? generator : zero(generator)
    observe_harmonics!(profile, dressed_generator[FE.triindex(n, 0)])

    for j in 1:n
      record_generator_node!(profile, K, dressed_generator, n, j)
      node = FE.dressed_generator_node(K, dressed_generator, n, j, generator)
      dressed_generator[FE.triindex(n, j)] = observe_harmonics!(profile, node)
    end

    for j in 1:n
      record_kick_node!(profile, K, Kdot, dressed_kick_derivative, n, j)
      node = FE.dressed_kick_derivative_node(
        K, Kdot, dressed_kick_derivative, n, j, generator
      )
      dressed_kick_derivative[FE.triindex(n, j)] = observe_harmonics!(profile, node)
    end

    record_assembly!(profile, dressed_generator, dressed_kick_derivative, n)
    raw_resolvent = FE.assemble_resolvent(
      dressed_generator, dressed_kick_derivative, n, generator
    )
    observe_harmonics!(profile, raw_resolvent)
    profile.periodic_simplify_calls += 1
    profile.simplify_harmonic_inputs += length(raw_resolvent)
    resolvent = observe_harmonics!(profile, FE.SQA.simplify(raw_resolvent))

    profile.component_simplify_calls += 1
    effective_n = FE.SQA.simplify(FE.time_average(resolvent))
    push!(effective, effective_n)

    if n < order - 1
      raw_kick = FE.antiderivative(FE.remove_average(resolvent), gauge)
      observe_harmonics!(profile, raw_kick)
      profile.periodic_simplify_calls += 1
      profile.simplify_harmonic_inputs += length(raw_kick)
      next_kick = observe_harmonics!(profile, FE.SQA.simplify(raw_kick))
      push!(K, next_kick)
      push!(Kdot, observe_harmonics!(profile, FE.derivative(next_kick)))
    end
  end

  return K, effective, profile
end

function verify_coherent_profile(H, ω, t, order)
  generator = FE.harmonics(FE.qadd(H), ω, t)
  K, effective, profile = coherent_work_profile(generator, VanVleck(), order)
  reference = floquet_expansion(generator, VanVleck(), order)
  K == getfield(reference, :kick_components) || error("profile kick mismatch")
  effective == getfield(reference, :effective_components) ||
    error("profile effective-component mismatch")
  return generator, profile
end

function print_coherent_profile(label, H, ω, t, order)
  generator, profile = verify_coherent_profile(H, ω, t, order)
  fields = (
    "workload=$(label)",
    "order=$(order)",
    "input_harmonics=$(length(generator))",
    "dressed_generator_nodes=$(profile.dressed_generator_nodes)",
    "dressed_kick_nodes=$(profile.dressed_kick_nodes)",
    "periodic_commutators=$(profile.periodic_commutators)",
    "harmonic_commutator_pairs=$(profile.harmonic_commutator_pairs)",
    "assembly_terms=$(profile.assembly_terms)",
    "assembly_harmonic_inputs=$(profile.assembly_harmonic_inputs)",
    "periodic_simplify_calls=$(profile.periodic_simplify_calls)",
    "component_simplify_calls=$(profile.component_simplify_calls)",
    "simplify_harmonic_inputs=$(profile.simplify_harmonic_inputs)",
    "materialized_harmonics=$(profile.materialized_harmonics)",
    "peak_harmonics=$(profile.peak_harmonics)",
  )
  println("COHERENT_PROFILE ", join(fields, " "))
  return nothing
end

function print_coherent_work_profiles()
  qubit_H, qubit_ω, qubit_t, _, _ = cp_hfe_qubit_workload()
  boson_H, boson_ω, boson_t, _ = cp_hfe_bosonic_workload()
  for order in 1:4
    print_coherent_profile("qubit", qubit_H, qubit_ω, qubit_t, order)
    print_coherent_profile("kerr", boson_H, boson_ω, boson_t, order)
  end
  return nothing
end

using BenchmarkTools: @benchmarkable
using Symbolics: @variables

const FE_BLOCH = FloquetExpansions

include(joinpath(@__DIR__, "..", "test", "helpers", "bloch_evaluation_plan.jl"))

bloch_operator_product(left, right) = left * right
bloch_hamiltonian_inverse(harmonic) = 1 // harmonic
bloch_simplify(value) = FE_BLOCH.SQA.simplify(value)

mutable struct BlochAlgebraProfile
  simplify_calls::Int
  simplify_term_inputs::Int
  simplify_term_outputs::Int
  peak_simplify_input_terms::Int
  peak_simplify_output_terms::Int
end

BlochAlgebraProfile() = BlochAlgebraProfile(0, 0, 0, 0, 0)

function profiled_bloch_simplify(profile::BlochAlgebraProfile, value)
  input_terms = length(value)
  profile.simplify_calls += 1
  profile.simplify_term_inputs += input_terms
  profile.peak_simplify_input_terms = max(profile.peak_simplify_input_terms, input_terms)

  simplified = bloch_simplify(value)
  output_terms = length(simplified)
  profile.simplify_term_outputs += output_terms
  profile.peak_simplify_output_terms = max(profile.peak_simplify_output_terms, output_terms)
  return simplified
end

function bloch_qubit_workload()
  pauli = PauliSpace(:bloch_benchmark_qubit)
  σx = Pauli(pauli, :sigma, 1)
  @variables ω_bloch_qubit::Real t_bloch_qubit::Real Ω_bloch_qubit::Real
  ω = ω_bloch_qubit
  t = t_bloch_qubit
  H = Ω_bloch_qubit * cos(ω * t) * σx
  return H, ω, t
end

function bloch_kerr_workload()
  fock = FockSpace(:bloch_benchmark_kerr)
  a = Destroy(fock, :a)
  @variables ω_bloch_kerr::Real t_bloch_kerr::Real Δ_bloch_kerr::Real
  @variables K_bloch_kerr::Real ε_bloch_kerr::Real
  ω = ω_bloch_kerr
  t = t_bloch_kerr
  number = a' * a
  H =
    Δ_bloch_kerr * number + K_bloch_kerr * a'^2 * a^2 + ε_bloch_kerr * cos(ω * t) * (a + a')
  return H, ω, t
end

function bloch_benchmark_context(H, ω, t, order)
  generator = FE_BLOCH.harmonics(FE_BLOCH.qadd(H), ω, t)
  support = collect(keys(generator))
  plan = compile_bloch_evaluation_plan(support, order)
  components = getfield(generator, :components)
  zero_component = getfield(generator, :zero_component)
  return generator, support, plan, components, zero_component
end

function evaluate_bloch_hamiltonian(plan, components, zero_component)
  return evaluate_bloch_evaluation_plan(
    plan,
    components;
    product=bloch_operator_product,
    inverse_weight=bloch_hamiltonian_inverse,
    zero_component,
    simplifier=bloch_simplify,
  )
end

function profile_bloch_hamiltonian(plan, components, zero_component)
  profile = BlochAlgebraProfile()
  simplifier = value -> profiled_bloch_simplify(profile, value)
  result = evaluate_bloch_evaluation_plan(
    plan,
    components;
    product=bloch_operator_product,
    inverse_weight=bloch_hamiltonian_inverse,
    zero_component,
    simplifier,
  )
  return result, profile
end

function solve_bloch_hamiltonian(H, ω, t, order)
  _, _, plan, components, zero_component = bloch_benchmark_context(H, ω, t, order)
  return evaluate_bloch_hamiltonian(plan, components, zero_component)
end

function bloch_input_term_count(generator)
  return sum(length(component) for component in values(getfield(generator, :components)))
end

function print_bloch_plan_profile(label, H, ω, t, order)
  generator, _, plan, components, zero_component = bloch_benchmark_context(H, ω, t, order)
  _, algebra = profile_bloch_hamiltonian(plan, components, zero_component)
  counts = plan.counts
  fields = (
    "workload=$(label)",
    "order=$(order)",
    "input_harmonics=$(length(generator))",
    "input_terms=$(bloch_input_term_count(generator))",
    "residual_nodes=$(counts.residual_nodes)",
    "wave_nodes=$(counts.wave_nodes)",
    "effective_nodes=$(counts.effective_nodes)",
    "generator_products=$(counts.generator_products)",
    "fold_products=$(counts.fold_products)",
    "total_products=$(counts.generator_products + counts.fold_products)",
    "simplify_calls=$(algebra.simplify_calls)",
    "simplify_term_inputs=$(algebra.simplify_term_inputs)",
    "simplify_term_outputs=$(algebra.simplify_term_outputs)",
    "peak_simplify_input_terms=$(algebra.peak_simplify_input_terms)",
    "peak_simplify_output_terms=$(algebra.peak_simplify_output_terms)",
  )
  println("BLOCH_PLAN_PROFILE ", join(fields, " "))
  return nothing
end

function print_bloch_plan_profiles()
  qubit_H, qubit_ω, qubit_t = bloch_qubit_workload()
  kerr_H, kerr_ω, kerr_t = bloch_kerr_workload()
  for order in 2:4
    print_bloch_plan_profile("qubit", qubit_H, qubit_ω, qubit_t, order)
    print_bloch_plan_profile("kerr", kerr_H, kerr_ω, kerr_t, order)
  end
  return nothing
end

function benchmark_bloch_evaluation!(suite)
  workloads = (
    "Driven qubit" => bloch_qubit_workload(),
    "Driven Kerr resonator" => bloch_kerr_workload(),
  )

  for (label, (H, ω, t)) in workloads, order in 2:4
    _, support, plan, components, zero_component = bloch_benchmark_context(H, ω, t, order)

    suite["Bloch Evaluation"][label]["order $order"]["plan build"] = @benchmarkable compile_bloch_evaluation_plan(
      $support, $order
    )
    suite["Bloch Evaluation"][label]["order $order"]["precompiled evaluation"] = @benchmarkable evaluate_bloch_hamiltonian(
      $plan, $components, $zero_component
    )
    suite["Bloch Evaluation"][label]["order $order"]["compile + evaluation"] = @benchmarkable solve_bloch_hamiltonian(
      $H, $ω, $t, $order
    )
    suite["Bloch Evaluation"][label]["order $order"]["production Hori-Deprit"] = @benchmarkable floquet_expansion(
      $H, $ω, $t, VanVleck(), $order
    )
  end
  return nothing
end

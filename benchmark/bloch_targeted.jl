using BenchmarkTools
using FloquetExpansions

include(joinpath(@__DIR__, "bloch_evaluation.jl"))

function print_trial(label, trial)
  estimate = BenchmarkTools.median(trial)
  println(
    "BLOCH_BENCHMARK ",
    label,
    " time_ns=", estimate.time,
    " memory_bytes=", estimate.memory,
    " allocs=", estimate.allocs,
  )
  return nothing
end

function run_targeted_bloch_benchmarks()
  workloads = (
    "qubit" => bloch_qubit_workload(),
    "kerr" => bloch_kerr_workload(),
  )

  print_bloch_plan_profiles()

  for (label, (H, ω, t)) in workloads, order in 2:4
    _, support, plan, components, zero_component = bloch_benchmark_context(H, ω, t, order)

    plan_trial = @benchmark compile_bloch_evaluation_plan($support, $order) seconds=0.4
    eval_trial = @benchmark evaluate_bloch_hamiltonian($plan, $components, $zero_component) seconds=0.4
    end_to_end_trial = @benchmark solve_bloch_hamiltonian($H, $ω, $t, $order) seconds=0.4
    hori_trial = @benchmark floquet_expansion($H, $ω, $t, VanVleck(), $order) seconds=0.4

    prefix = "workload=$(label) order=$(order)"
    print_trial("$(prefix) stage=plan_build", plan_trial)
    print_trial("$(prefix) stage=precompiled_evaluation", eval_trial)
    print_trial("$(prefix) stage=compile_evaluation", end_to_end_trial)
    print_trial("$(prefix) stage=hori_deprit", hori_trial)
  end
  return nothing
end

run_targeted_bloch_benchmarks()

using BenchmarkTools: @benchmarkable
using Symbolics: @variables

include(joinpath(@__DIR__, "..", "test", "helpers", "lyndon_log_evaluator.jl"))

const FE_LLB = FloquetExpansions

llb_product(left, right) = left * right
llb_simplify(value) = FE_LLB.SQA.simplify(value)

function llb_qubit_workload()
  pauli = PauliSpace(:lyndon_benchmark_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω_llb_q::Real t_llb_q::Real Ω_llb_q::Real Δ_llb_q::Real
  H = Δ_llb_q * σz + Ω_llb_q * cos(ω_llb_q * t_llb_q) * σx
  return H, ω_llb_q, t_llb_q
end

function llb_kerr_workload()
  fock = FockSpace(:lyndon_benchmark_kerr)
  a = Destroy(fock, :a)
  @variables ω_llb_k::Real t_llb_k::Real Δ_llb_k::Real K_llb_k::Real ε_llb_k::Real
  number = a' * a
  H = Δ_llb_k * number + K_llb_k * a'^2 * a^2 + ε_llb_k * cos(ω_llb_k * t_llb_k) * (a + a')
  return H, ω_llb_k, t_llb_k
end

function llb_mercator_replay(
  normalized::Vector{Dict{H,T}}; product=llb_product, simplifier=llb_simplify
) where {H,T}
  order = length(normalized)
  counts = FE_LLB.BlochVanVleckCounts()
  powers = [[Dict{H,T}() for _ in 1:max(order, 1)] for _ in 1:max(order, 1)]
  log_embedding = Vector{Dict{H,T}}()

  for n in 1:order
    powers[1][n] = normalized[n]
    generator_n = copy(normalized[n])
    for power in 2:n
      coefficient = Dict{H,T}()
      for k in 1:(n - power + 1)
        counts.log_products += 1
        contribution = FE_LLB.bloch_vv_periodic_product(
          normalized[k], powers[power - 1][n - k], product, simplifier, counts
        )
        coefficient = FE_LLB.bloch_vv_add(coefficient, contribution, simplifier)
      end
      powers[power][n] = coefficient
      weight = (-1)^(power + 1) * (1 // power)
      generator_n = FE_LLB.bloch_vv_add(
        generator_n, FE_LLB.bloch_vv_scale(weight, coefficient, simplifier), simplifier
      )
    end
    push!(log_embedding, generator_n)
  end
  return log_embedding, counts
end

function llb_context(H, ω, t, order)
  generator = FE_LLB.harmonics(FE_LLB.qadd(H), ω, t)
  components = Dict(harmonic => generator[harmonic] for harmonic in keys(generator))
  zero_component = generator.zero_component
  bloch_plan = FE_LLB.compile_bloch_projection_plan(keys(components), order)
  bloch = FE_LLB.evaluate_bloch_projection_plan(
    bloch_plan,
    components;
    product=llb_product,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component,
    simplifier=llb_simplify,
  )
  direct = FE_LLB.bloch_van_vleck_reconstruction(
    bloch_plan, bloch; product=llb_product, zero_component, simplifier=llb_simplify
  )
  lyndon = compile_lyndon_log_evaluation_plan(keys(components), order)
  _, mercator_counts = llb_mercator_replay(direct.normalized_embedding)
  return components, zero_component, direct, lyndon, mercator_counts
end

function print_lyndon_backend_profile(label, H, ω, t, order)
  _, _, direct, lyndon, mercator = llb_context(H, ω, t, order)
  println(
    "LYNDON_BACKEND_PROFILE ",
    "workload=$(label) ",
    "order=$(order) ",
    "mercator_series_products=$(mercator.log_products) ",
    "mercator_backend_products=$(mercator.harmonic_products) ",
    "lyndon_bracket_nodes=$(length(lyndon.brackets)) ",
    "lyndon_backend_products=$(lyndon_backend_products(lyndon)) ",
    "canonical_log_harmonics=$(sum(length, direct.log_embedding))",
  )
  return nothing
end

function benchmark_lyndon_log_evaluation!(suite)
  workloads = (
    "Driven qubit" => llb_qubit_workload(), "Driven Kerr resonator" => llb_kerr_workload()
  )

  for (label, (H, ω, t)) in workloads, order in 2:4
    components, zero_component, direct, lyndon, _ = llb_context(H, ω, t, order)
    support = collect(keys(components))
    normalized = direct.normalized_embedding

    suite["Connected Log Evaluation"][label]["order $order"]["Mercator replay"] = @benchmarkable llb_mercator_replay(
      $normalized
    )
    suite["Connected Log Evaluation"][label]["order $order"]["Lyndon DAG"] = @benchmarkable evaluate_lyndon_log_plan(
      $lyndon,
      $components;
      product=llb_product,
      zero_component=($zero_component),
      simplifier=llb_simplify,
    )
    suite["Connected Log Evaluation"][label]["order $order"]["Lyndon plan compile"] = @benchmarkable compile_lyndon_log_evaluation_plan(
      $support, $order
    )

    print_lyndon_backend_profile(label, H, ω, t, order)
  end
  return nothing
end

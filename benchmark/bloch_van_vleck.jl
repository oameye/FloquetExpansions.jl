using BenchmarkTools: @benchmarkable
using Symbolics: @variables

const FE_BVVB = FloquetExpansions

include(joinpath(@__DIR__, "..", "test", "helpers", "bloch_van_vleck_projection.jl"))

bvv_product(left, right) = left * right
bvv_inverse_weight(harmonic) = 1 // harmonic
bvv_simplify(value) = FE_BVVB.SQA.simplify(value)

function bvv_qubit_workload()
  pauli = PauliSpace(:bvv_benchmark_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω_bvv_qubit::Real t_bvv_qubit::Real Ω_bvv_qubit::Real Δ_bvv_qubit::Real
  H = Δ_bvv_qubit * σz + Ω_bvv_qubit * cos(ω_bvv_qubit * t_bvv_qubit) * σx
  return H, ω_bvv_qubit, t_bvv_qubit
end

function bvv_kerr_workload()
  fock = FockSpace(:bvv_benchmark_kerr)
  a = Destroy(fock, :a)
  @variables ω_bvv_kerr::Real t_bvv_kerr::Real Δ_bvv_kerr::Real
  @variables K_bvv_kerr::Real ε_bvv_kerr::Real
  number = a' * a
  H =
    Δ_bvv_kerr * number +
    K_bvv_kerr * a'^2 * a^2 +
    ε_bvv_kerr * cos(ω_bvv_kerr * t_bvv_kerr) * (a + a')
  return H, ω_bvv_kerr, t_bvv_kerr
end

function bvv_generator(H, ω, t)
  return FE_BVVB.harmonics(FE_BVVB.qadd(H), ω, t)
end

function bvv_projection(generator, order)
  return evaluate_periodic_bloch_projection(
    generator,
    order;
    product=bvv_product,
    inverse_weight=bvv_inverse_weight,
    simplifier=bvv_simplify,
  )
end

function bvv_reconstruction(bloch, generator)
  return bloch_van_vleck_projection(
    bloch, generator; product=bvv_product, simplifier=bvv_simplify
  )
end

function bvv_end_to_end(H, ω, t, order)
  generator = bvv_generator(H, ω, t)
  bloch = bvv_projection(generator, order)
  return bvv_reconstruction(bloch, generator)
end

function bvv_benchmark_context(H, ω, t, order)
  generator = bvv_generator(H, ω, t)
  plan = FE_BVVB.compile_bloch_projection_plan(collect(keys(generator)), order)
  bloch = bvv_projection(generator, order)
  converted = bvv_reconstruction(bloch, generator)
  return generator, plan, bloch, converted
end

function print_bvv_profile(label, H, ω, t, order)
  _, plan, _, converted = bvv_benchmark_context(H, ω, t, order)
  projection = plan.counts
  conversion = converted.counts
  fields = (
    "workload=$(label)",
    "order=$(order)",
    "projection_generator_products=$(projection.generator_products)",
    "projection_fold_products=$(projection.fold_products)",
    "factor_products=$(conversion.factor_products)",
    "log_products=$(conversion.log_products)",
    "expected_log_products=$(expected_projection_mercator_products(order - 1))",
    "inverse_products=$(conversion.inverse_products)",
    "similarity_products=$(conversion.similarity_products)",
    "conversion_component_products=$(conversion.periodic_component_products)",
  )
  println("BLOCH_VV_PROFILE ", join(fields, " "))
  return nothing
end

function print_bvv_profiles()
  qubit_H, qubit_ω, qubit_t = bvv_qubit_workload()
  kerr_H, kerr_ω, kerr_t = bvv_kerr_workload()
  for order in 2:4
    print_bvv_profile("qubit", qubit_H, qubit_ω, qubit_t, order)
    print_bvv_profile("kerr", kerr_H, kerr_ω, kerr_t, order)
  end
  return nothing
end

function benchmark_bloch_van_vleck!(suite)
  workloads = (
    "Driven qubit" => bvv_qubit_workload(),
    "Driven Kerr resonator" => bvv_kerr_workload(),
  )

  for (label, (H, ω, t)) in workloads, order in 2:4
    generator, _, bloch, converted = bvv_benchmark_context(H, ω, t, order)
    normalized = converted.normalized_embedding
    static_factor = converted.static_factor
    effective = bloch.effective

    suite["Bloch Canonical Reconstruction"][label]["order $order"]["projection"] =
      @benchmarkable bvv_projection($generator, $order)
    suite["Bloch Canonical Reconstruction"][label]["order $order"]["canonical reconstruction"] =
      @benchmarkable bvv_reconstruction($bloch, $generator)
    suite["Bloch Canonical Reconstruction"][label]["order $order"]["Mercator log replay"] =
      @benchmarkable replay_mercator_log(
        $normalized, $generator; product=bvv_product, simplifier=bvv_simplify
      )
    suite["Bloch Canonical Reconstruction"][label]["order $order"]["static similarity replay"] =
      @benchmarkable replay_static_similarity(
        $effective, $static_factor; product=bvv_product, simplifier=bvv_simplify
      )
    suite["Bloch Canonical Reconstruction"][label]["order $order"]["end to end"] =
      @benchmarkable bvv_end_to_end($H, $ω, $t, $order)
    suite["Bloch Canonical Reconstruction"][label]["order $order"]["production Hori-Deprit"] =
      @benchmarkable floquet_expansion($H, $ω, $t, VanVleck(), $order)
  end
  return nothing
end

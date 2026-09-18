using BenchmarkTools: @benchmarkable

include(
  joinpath(@__DIR__, "..", "test", "helpers", "connected_van_vleck_reconstruction.jl")
)

const FE_CVVB = FloquetExpansions

cvvb_product(left, right) = left * right
cvvb_simplify(value) = FE_CVVB.SQA.simplify(value)

function cvvb_context(H, ω, t, order)
  generator = FE_CVVB.harmonics(FE_CVVB.qadd(H), ω, t)
  components = Dict(harmonic => generator[harmonic] for harmonic in keys(generator))
  zero_component = generator.zero_component
  projection_plan = FE_CVVB.compile_bloch_projection_plan(keys(components), order)
  bloch = FE_CVVB.evaluate_bloch_projection_plan(
    projection_plan,
    components;
    product=cvvb_product,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component,
    simplifier=cvvb_simplify,
  )
  direct = FE_CVVB.bloch_van_vleck_reconstruction(
    projection_plan, bloch; product=cvvb_product, zero_component, simplifier=cvvb_simplify
  )
  log_plan = compile_lyndon_log_evaluation_plan(keys(components), order)
  static_plan = compile_static_sector_exp_plan(log_plan, projection_plan.zero_harmonic)
  connected = connected_van_vleck_reconstruction(
    projection_plan,
    bloch,
    components,
    log_plan,
    static_plan;
    product=cvvb_product,
    zero_component,
    simplifier=cvvb_simplify,
  )
  return (
    components,
    zero_component,
    projection_plan,
    bloch,
    direct,
    log_plan,
    static_plan,
    connected,
  )
end

function cvvb_direct_reconstruction(projection_plan, bloch, zero_component)
  return FE_CVVB.bloch_van_vleck_reconstruction(
    projection_plan, bloch; product=cvvb_product, zero_component, simplifier=cvvb_simplify
  )
end

function cvvb_connected_reconstruction(
  projection_plan, bloch, components, log_plan, static_plan, zero_component
)
  return connected_van_vleck_reconstruction(
    projection_plan,
    bloch,
    components,
    log_plan,
    static_plan;
    product=cvvb_product,
    zero_component,
    simplifier=cvvb_simplify,
  )
end

function cvvb_direct_core(projection_plan, components, zero_component)
  bloch = FE_CVVB.evaluate_bloch_projection_plan(
    projection_plan,
    components;
    product=cvvb_product,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component,
    simplifier=cvvb_simplify,
  )
  return cvvb_direct_reconstruction(projection_plan, bloch, zero_component)
end

function cvvb_connected_core(
  projection_plan, log_plan, static_plan, components, zero_component
)
  bloch = FE_CVVB.evaluate_bloch_projection_plan(
    projection_plan,
    components;
    product=cvvb_product,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component,
    simplifier=cvvb_simplify,
  )
  return cvvb_connected_reconstruction(
    projection_plan, bloch, components, log_plan, static_plan, zero_component
  )
end

function print_connected_reconstruction_profile(label, H, ω, t, order)
  _, _, _, _, direct, log_plan, static_plan, connected = cvvb_context(H, ω, t, order)
  direct_products =
    direct.counts.harmonic_products +
    direct.counts.inverse_products +
    direct.counts.similarity_products
  connected_products = connected_reconstruction_backend_products(connected, log_plan)
  println(
    "CONNECTED_RECONSTRUCTION_PROFILE ",
    "workload=$(label) ",
    "order=$(order) ",
    "direct_backend_products=$(direct_products) ",
    "connected_backend_products=$(connected_products) ",
    "connected_bracket_products=$(lyndon_backend_products(log_plan)) ",
    "connected_exp_harmonic_products=$(connected.counts.harmonic_products) ",
    "compiled_exp_harmonic_products=$(static_plan.product_count) ",
    "connected_static_products=$(connected.counts.inverse_products + connected.counts.similarity_products)",
  )
  return nothing
end

function benchmark_connected_van_vleck_reconstruction!(suite)
  workloads = (
    "Driven qubit" => llb_qubit_workload(), "Driven Kerr resonator" => llb_kerr_workload()
  )

  for (label, (H, ω, t)) in workloads, order in 2:4
    (
      components,
      zero_component,
      projection_plan,
      bloch,
      _,
      log_plan,
      static_plan,
      connected,
    ) = cvvb_context(H, ω, t, order)
    log_embedding = connected.log_embedding
    identity_component = one(first(bloch.effective))
    zero_harmonic = projection_plan.zero_harmonic

    suite["Connected Canonical Reconstruction"][label]["order $order"]["direct reconstruction"] = @benchmarkable cvvb_direct_reconstruction(
      $projection_plan, $bloch, $zero_component
    )
    suite["Connected Canonical Reconstruction"][label]["order $order"]["connected reconstruction"] = @benchmarkable cvvb_connected_reconstruction(
      $projection_plan,
      $bloch,
      $components,
      $log_plan,
      $static_plan,
      $zero_component,
    )
    suite["Connected Canonical Reconstruction"][label]["order $order"]["static full periodic exp"] = @benchmarkable connected_static_factor(
      $log_embedding,
      $zero_harmonic,
      $identity_component,
      $zero_component;
      product=cvvb_product,
      simplifier=cvvb_simplify,
    )
    suite["Connected Canonical Reconstruction"][label]["order $order"]["static pruned exp"] = @benchmarkable evaluate_static_sector_exp_plan(
      $static_plan,
      $log_embedding,
      $identity_component,
      $zero_component;
      product=cvvb_product,
      simplifier=cvvb_simplify,
    )
    suite["Connected Canonical Reconstruction"][label]["order $order"]["direct core"] = @benchmarkable cvvb_direct_core(
      $projection_plan, $components, $zero_component
    )
    suite["Connected Canonical Reconstruction"][label]["order $order"]["connected core"] = @benchmarkable cvvb_connected_core(
      $projection_plan, $log_plan, $static_plan, $components, $zero_component
    )

    print_connected_reconstruction_profile(label, H, ω, t, order)
  end
  return nothing
end

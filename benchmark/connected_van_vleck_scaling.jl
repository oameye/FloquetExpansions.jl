using FloquetExpansions
using Statistics: median

include("workflows.jl")
include("lyndon_log_evaluation.jl")
include("connected_van_vleck_reconstruction.jl")

length(ARGS) == 2 || error("usage: connected_van_vleck_scaling.jl <qubit|kerr> <order>")
workload_name = ARGS[1]
order = parse(Int, ARGS[2])
order >= 4 || error("scaling benchmark is intended for order >= 4")

H, ω, t = if workload_name == "qubit"
  llb_qubit_workload()
elseif workload_name == "kerr"
  llb_kerr_workload()
else
  error("unknown workload: $(workload_name)")
end

function scaling_measure(f; samples=3)
  f()
  times = Float64[]
  for _ in 1:samples
    GC.gc()
    push!(times, @elapsed f())
  end
  GC.gc()
  memory = @allocated f()
  return median(times), memory
end

(
  components,
  zero_component,
  projection_plan,
  bloch,
  direct,
  log_plan,
  static_plan,
  connected,
) = cvvb_context(H, ω, t, order)

print_connected_reconstruction_profile(workload_name, H, ω, t, order)
println(
  "CONNECTED_SCALING_STRUCTURE ",
  "workload=$(workload_name) ",
  "order=$(order) ",
  "lyndon_brackets=$(length(log_plan.brackets)) ",
  "lyndon_backend_products=$(lyndon_backend_products(log_plan)) ",
  "static_exp_backend_products=$(static_plan.product_count) ",
  "connected_backend_products=$(connected_reconstruction_backend_products(connected, log_plan))",
)

direct_reconstruction =
  () -> cvvb_direct_reconstruction(projection_plan, bloch, zero_component)
connected_reconstruction =
  () -> cvvb_connected_reconstruction(
    projection_plan, bloch, components, log_plan, static_plan, zero_component
  )
direct_core = () -> cvvb_direct_core(projection_plan, components, zero_component)
connected_core =
  () -> cvvb_connected_core(
    projection_plan, log_plan, static_plan, components, zero_component
  )
lyndon_plan_compile = () -> compile_lyndon_log_evaluation_plan(keys(components), order)
static_plan_compile =
  () -> compile_static_sector_exp_plan(log_plan, projection_plan.zero_harmonic)
production = () -> begin
  expansion = floquet_expansion(H, ω, t, VanVleck(), order)
  effective_generator(expansion), micromotion(expansion)
end

for (label, f) in (
  "direct_reconstruction" => direct_reconstruction,
  "connected_reconstruction" => connected_reconstruction,
  "direct_core" => direct_core,
  "connected_core" => connected_core,
  "lyndon_plan_compile" => lyndon_plan_compile,
  "static_plan_compile" => static_plan_compile,
  "production_hori_deprit" => production,
)
  seconds, memory = scaling_measure(f)
  println(
    "CONNECTED_SCALING_RESULT ",
    "workload=$(workload_name) ",
    "order=$(order) ",
    "stage=$(label) ",
    "seconds=$(seconds) ",
    "memory=$(memory)",
  )
end

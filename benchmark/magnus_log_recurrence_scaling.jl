using FloquetExpansions
using Statistics: median

include("lyndon_log_evaluation.jl")
include(joinpath(@__DIR__, "..", "test", "helpers", "magnus_log_recurrence.jl"))

length(ARGS) == 2 || error("usage: magnus_log_recurrence_scaling.jl <qubit|kerr> <order>")
workload_name = ARGS[1]
order = parse(Int, ARGS[2])

H, ω, t = if workload_name == "qubit"
  llb_qubit_workload()
elseif workload_name == "kerr"
  llb_kerr_workload()
else
  error("unknown workload: $(workload_name)")
end

function mlr_measure(f; samples=3)
  f()
  times = Float64[]
  for _ in 1:samples
    GC.gc()
    push!(times, @elapsed f())
  end
  GC.gc()
  return median(times), @allocated(f())
end

components, zero_component, direct, lyndon, mercator_counts = llb_context(H, ω, t, order)
normalized = direct.normalized_embedding
magnus, magnus_counts = connected_log_magnus_recurrence(
  normalized; product=llb_product, zero_component, simplifier=llb_simplify
)
@assert all(
  begin
    harmonics = union(keys(magnus[n]), keys(direct.log_embedding[n]))
    all(
      iszero(
        llb_simplify(
          get(magnus[n], harmonic, zero_component) -
          get(direct.log_embedding[n], harmonic, zero_component)
        )
      ) for harmonic in harmonics
    )
  end for n in eachindex(magnus)
)

println(
  "MAGNUS_LOG_STRUCTURE ",
  "workload=$(workload_name) ",
  "order=$(order) ",
  "mercator_products=$(mercator_counts.harmonic_products) ",
  "lyndon_products=$(lyndon_backend_products(lyndon)) ",
  "magnus_products=$(magnus_counts.harmonic_products)",
)

for (label, f) in (
  "mercator" => () -> llb_mercator_replay(normalized),
  "lyndon_eval" => () -> evaluate_lyndon_log_plan(
    lyndon,
    components;
    product=llb_product,
    zero_component,
    simplifier=llb_simplify,
  ),
  "lyndon_compile" => () -> compile_lyndon_log_evaluation_plan(keys(components), order),
  "magnus_recurrence" => () -> connected_log_magnus_recurrence(
    normalized; product=llb_product, zero_component, simplifier=llb_simplify
  ),
)
  seconds, memory = mlr_measure(f)
  println(
    "MAGNUS_LOG_RESULT ",
    "workload=$(workload_name) ",
    "order=$(order) ",
    "stage=$(label) ",
    "seconds=$(seconds) ",
    "memory=$(memory)",
  )
end

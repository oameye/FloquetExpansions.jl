using FloquetExpansions
using Statistics: median

include("workflows.jl")
include(joinpath(@__DIR__, "..", "test", "helpers", "magnus_van_vleck_reconstruction.jl"))

const FE_MVVB = FloquetExpansions

length(ARGS) == 2 || error("usage: magnus_van_vleck_scaling.jl <qubit|kerr> <order>")
workload_name = ARGS[1]
order = parse(Int, ARGS[2])

H, ω, t = if workload_name == "qubit"
  llb_qubit_workload()
elseif workload_name == "kerr"
  llb_kerr_workload()
else
  error("unknown workload: $(workload_name)")
end

mvvb_product(left, right) = left * right
mvvb_simplify(value) = FE_MVVB.SQA.simplify(value)

function mvvb_measure(f; samples=3)
  f()
  times = Float64[]
  for _ in 1:samples
    GC.gc()
    push!(times, @elapsed f())
  end
  GC.gc()
  return median(times), @allocated(f())
end

generator = FE_MVVB.harmonics(FE_MVVB.qadd(H), ω, t)
components = Dict(harmonic => generator[harmonic] for harmonic in keys(generator))
zero_component = generator.zero_component
projection_plan = FE_MVVB.compile_bloch_projection_plan(keys(components), order)
bloch = FE_MVVB.evaluate_bloch_projection_plan(
  projection_plan,
  components;
  product=mvvb_product,
  inverse_weight=harmonic -> 1 // harmonic,
  zero_component,
  simplifier=mvvb_simplify,
)
direct = FE_MVVB.bloch_van_vleck_reconstruction(
  projection_plan, bloch; product=mvvb_product, zero_component, simplifier=mvvb_simplify
)
magnus = magnus_van_vleck_reconstruction(
  projection_plan, bloch; product=mvvb_product, zero_component, simplifier=mvvb_simplify
)
@assert all(
  iszero(mvvb_simplify(magnus.effective[n] - direct.effective[n])) for
  n in eachindex(direct.effective)
)

println(
  "MAGNUS_VV_STRUCTURE ",
  "workload=$(workload_name) ",
  "order=$(order) ",
  "direct_harmonic_products=$(direct.counts.harmonic_products) ",
  "magnus_harmonic_products=$(magnus.counts.harmonic_products)",
)

direct_reconstruction = () -> FE_MVVB.bloch_van_vleck_reconstruction(
  projection_plan, bloch; product=mvvb_product, zero_component, simplifier=mvvb_simplify
)
magnus_reconstruction = () -> magnus_van_vleck_reconstruction(
  projection_plan, bloch; product=mvvb_product, zero_component, simplifier=mvvb_simplify
)
direct_core = () -> begin
  local_bloch = FE_MVVB.evaluate_bloch_projection_plan(
    projection_plan,
    components;
    product=mvvb_product,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component,
    simplifier=mvvb_simplify,
  )
  FE_MVVB.bloch_van_vleck_reconstruction(
    projection_plan,
    local_bloch;
    product=mvvb_product,
    zero_component,
    simplifier=mvvb_simplify,
  )
end
magnus_core = () -> begin
  local_bloch = FE_MVVB.evaluate_bloch_projection_plan(
    projection_plan,
    components;
    product=mvvb_product,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component,
    simplifier=mvvb_simplify,
  )
  magnus_van_vleck_reconstruction(
    projection_plan,
    local_bloch;
    product=mvvb_product,
    zero_component,
    simplifier=mvvb_simplify,
  )
end
production = () -> begin
  expansion = floquet_expansion(H, ω, t, VanVleck(), order)
  effective_generator(expansion), micromotion(expansion)
end

for (label, f) in (
  "direct_reconstruction" => direct_reconstruction,
  "magnus_reconstruction" => magnus_reconstruction,
  "direct_core" => direct_core,
  "magnus_core" => magnus_core,
  "production_hori_deprit" => production,
)
  seconds, memory = mvvb_measure(f)
  println(
    "MAGNUS_VV_RESULT ",
    "workload=$(workload_name) ",
    "order=$(order) ",
    "stage=$(label) ",
    "seconds=$(seconds) ",
    "memory=$(memory)",
  )
end

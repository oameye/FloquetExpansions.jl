using Statistics: median

include(joinpath(@__DIR__, "..", "test", "helpers", "lyndon_log_evaluator.jl"))

length(ARGS) == 1 || error("usage: lyndon_plan_compile_profile.jl <order>")
order = parse(Int, only(ARGS))
support = [-1, 0, 1]

function measure_stage(f; samples=5)
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

_, _, converted = connected_word_reconstruction(support, order)

function decomposition_stage()
  result = Vector{Dict{Int,Dict{Tuple,Rational{Int}}}}()
  for embedding in converted.log_embedding
    compiled = Dict{Int,Dict{Tuple,Rational{Int}}}()
    for (harmonic, polynomial) in embedding
      compiled[harmonic] = lyndon_decomposition(polynomial)
    end
    push!(result, compiled)
  end
  return result
end

for (label, f) in (
  "word_reconstruction" => (() -> connected_word_reconstruction(support, order)),
  "lyndon_decomposition" => decomposition_stage,
  "full_plan_compile" => (() -> compile_lyndon_log_evaluation_plan(support, order)),
)
  seconds, memory = measure_stage(f)
  println(
    "LYNDON_COMPILE_PROFILE ",
    "order=$(order) ",
    "stage=$(label) ",
    "seconds=$(seconds) ",
    "memory=$(memory)",
  )
end

plan = compile_lyndon_log_evaluation_plan(support, order)
println(
  "LYNDON_COMPILE_STRUCTURE ",
  "order=$(order) ",
  "brackets=$(length(plan.brackets)) ",
  "output_terms=$(sum(length(terms) for embedding in plan.outputs for terms in values(embedding)))",
)

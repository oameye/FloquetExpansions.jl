using BenchmarkTools
using FloquetExpansions

const SUITE = BenchmarkGroup()

include("kerr_resonator.jl")

benchmark_kerr_resonator!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))

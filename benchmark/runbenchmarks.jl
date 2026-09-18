using BenchmarkTools
using FloquetExpansions

const SUITE = BenchmarkGroup()

include("workflows.jl")
include("bloch_van_vleck.jl")

benchmark_fourier_expansion!(SUITE)
benchmark_floquet_expansion!(SUITE)
benchmark_positive_completion!(SUITE)
benchmark_bloch_van_vleck!(SUITE)
print_bvv_profiles()

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)
display(median(results))

BenchmarkTools.save(joinpath(@__DIR__, "benchmarks_output.json"), median(results))

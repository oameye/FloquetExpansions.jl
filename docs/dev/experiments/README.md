# Experiments

Throwaway scripts that established facts recorded in `../RESEARCH.md`. Not part of the package and
not run by the test suite. Kept because each one cost real effort and each is worth re-running
before trusting the claim it supports.

They need `SecondQuantizedAlgebra`, and additionally:

    collector_prototype.jl   SymbolicUtils, Symbolics
    gauge_experiment.jl      LinearAlgebra, Printf   (stdlib)
    nikolaev_experiment.jl   LinearAlgebra, Printf   (stdlib)

Run in a scratch environment, not the package environment:

    julia --project=/tmp/vv-scratch -e 'using Pkg
        Pkg.add(url="https://github.com/qojulia/SecondQuantizedAlgebra.jl", rev="main")
        Pkg.add(["SymbolicUtils","Symbolics"])'
    julia --project=/tmp/vv-scratch collector_prototype.jl

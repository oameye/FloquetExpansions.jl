# `PeriodicGenerator` is the Fourier boundary

`PeriodicGenerator` is the canonical representation of a periodic generator: it carries the drive frequency and its nonzero integer harmonics, while missing harmonics evaluate as zero. It replaces the Hamiltonian-only `PeriodicOperator`; no alias or separate periodic Liouvillian container is maintained.

## Gate

`make test` runs the testsets that hold this decision:

- `test/periodic_operator.jl`: "the drive frequency is part of the periodic generator"
- `test/periodic_operator.jl`: "constructors normalize zero harmonics"
- `test/liouvillian.jl`: "a zero periodic Liouvillian keeps its public zero prototype"

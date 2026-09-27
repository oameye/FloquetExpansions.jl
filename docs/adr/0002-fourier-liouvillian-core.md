# One expansion core for all generator algebras

Hamiltonian and Liouvillian inputs use the same `PeriodicGenerator{T}` Fourier representation and van Vleck recursion, so one `FloquetExpansion` model serves both. The component type carries the algebra through addition, commutators, differentiation, and simplification; no parallel Liouvillian expansion is maintained.

## Gate

`make test` runs the testsets that hold this decision:

- `test/generators/liouvillian.jl`: "Liouvillians use the common van Vleck expansion"
- `test/expansion/hori_deprit.jl`: "the parsing entry point agrees with the harmonic one"

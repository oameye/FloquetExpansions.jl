# The van Vleck gauge carries a public expansion-algorithm selector

`VanVleck` is parameterized by an `ExpansionAlgorithm`, chosen as `VanVleck(; algorithm=...)` with `HoriDeprit()` as the default and `BlochFeshbach()` as the alternative. Both algorithms compute the same retained effective components and micromotion, so the selector changes cost and intermediate quantities, never the result. The selector is public because neither algorithm dominates. The evidence is steady-state runtime, the BenchmarkTools median of `floquet_expansion` on a prepared `PeriodicGenerator`, measured on 2026-09-26 with Julia 1.13.0 on the workloads of `benchmark/workflows.jl` plus a Pauli qubit:

- On driven-qubit Hamiltonians the faster algorithm depends on the order and the harmonic content. At order 6 either wins, by at most a factor of 2. At order 10 Hori–Deprit is 11 times faster on the `NLevelSpace` qubit and 56 times faster on the Pauli qubit with harmonics ``0, \pm 1``.
- On the driven Kerr oscillator at order 6, Bloch/Feshbach is about 5 times faster.
- On the driven dissipative qubit at orders 4 to 6, Bloch/Feshbach is 14 to 21 times faster.

On Hamiltonians at high order, compiling the word-level plan takes 93 % of the Bloch/Feshbach time at order 8 and more than 99 % at order 10. The plan depends only on the harmonic support and order, so the Hamiltonian crossover would move if the plan were reused across calls.

## Considered options

- **Internal Bloch/Feshbach only**, as #97 proposed: no public selector, with Bloch/Feshbach kept as a tested internal backend. Rejected because the measured crossover depends on the generator, so no single internal choice is best.
- **Automatic selection by generator type**: Bloch/Feshbach for Liouvillians, Hori–Deprit for Hamiltonians. Rejected for now because the evidence is a handful of workloads, too thin to fix a rule that users cannot override.
- **An `algorithm` keyword on `floquet_expansion`** instead of a gauge parameter. Rejected because the gauge is already the value a `FloquetExpansion` stores about how it was produced, so a gauge parameter records the algorithm on the result without a second field or a second dispatch argument through `floquet_expansion`.

## Consequences

The algorithm is part of the gauge's type, so `VanVleck()` and `VanVleck(; algorithm=BlochFeshbach())` compare unequal and yield different `FloquetExpansion` types. Algorithm equivalence is a statement about retained coefficients, not about result types.

Algorithm-specific intermediate data, such as the wave operator and the Bloch effective generator, are not retained on the result. This differs from positive completion, which keeps its data behind `factorization` (ADR 0009), and holds until a consumer of expansion-algorithm diagnostics exists.

## Gate

`make test` runs the testsets in `test/expansion/expansion_algorithms.jl` that hold this decision:

- "Van Vleck algorithm selectors"
- "Bloch Feshbach matches Hamiltonian Hori Deprit Van Vleck"
- "Bloch Feshbach matches Liouvillian Hori Deprit"

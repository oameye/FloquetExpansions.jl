# The van Vleck gauge carries a public expansion-algorithm selector

`VanVleck` is parameterized by an `ExpansionAlgorithm`, chosen as `VanVleck(; algorithm=...)` with `HoriDeprit()` as the default and `BlochFeshbach()` as the alternative. Both algorithms compute the same retained effective components and micromotion, so the selector changes cost and intermediate quantities, never the result. The selector is public because neither algorithm dominates: on the benchmarks behind this decision, Hori–Deprit was 2 to 25 times faster on driven-qubit Hamiltonians at orders 6 to 10, while Bloch/Feshbach was about 12 times faster on a driven Kerr oscillator at order 6 and 11 to 20 times faster on qubit and cavity Liouvillians at orders 4 to 7.

## Considered options

- **Internal Bloch/Feshbach only**, as #97 proposed: no public selector, with Bloch/Feshbach kept as a tested internal backend. Rejected because the measured crossover depends on the generator, so no single internal choice is best.
- **Automatic selection by generator type**: Bloch/Feshbach for Liouvillians, Hori–Deprit for Hamiltonians. Rejected for now because the evidence is a handful of workloads, too thin to fix a rule that users cannot override.
- **An `algorithm` keyword on `floquet_expansion`** instead of a gauge parameter. Rejected because the gauge is already the value a `FloquetExpansion` stores about how it was produced, and every accessor and completion method already dispatches on `VanVleck`.

## Consequences

The algorithm is part of the gauge's type, so `VanVleck()` and `VanVleck(; algorithm=BlochFeshbach())` compare unequal and yield different `FloquetExpansion` types. Algorithm equivalence is a statement about retained coefficients, not about result types.

Algorithm-specific intermediate data, such as the wave operator and the Bloch effective generator, are not retained on the result. This differs from positive completion, which keeps its data behind `factorization` (ADR 0009), and holds until a consumer of expansion-algorithm diagnostics exists.

## Gate

`make test` runs the testsets in `test/expansion_algorithms.jl` that hold this decision:

- "Van Vleck algorithm selectors"
- "Bloch Feshbach matches Hamiltonian Hori Deprit Van Vleck"
- "Bloch Feshbach matches Liouvillian Hori Deprit"

# Architecture Map

Read this guide when a change crosses module boundaries, changes a symbolic representation, or changes the public expansion model. Read the relevant ADRs for the rationale and constraints behind the map.

## Data flow

```text
time-dependent SQA expression
        │
        ▼
PeriodicGenerator{T}  ──►  Van Vleck recursion  ──►  FloquetExpansion
        │                                             │
        └──────────────► QuasienergyOperator         ├─ effective_generator
                                                      └─ micromotion
```

`PeriodicGenerator{T}` is the shared Fourier boundary. Its component type carries the algebra through addition, commutators, derivatives, antiderivatives, and simplification. The expansion engine should remain independent of whether `T` is `SQA.QAdd` or `Liouvillian`.

## Module ownership

| Module | Owns | Public seam |
| --- | --- | --- |
| `periodic_operator.jl` | Fourier harmonics, drive frequency, harmonic calculus, gauges | `PeriodicGenerator`, `harmonics`, `time_average`, `derivative`, `antiderivative` |
| `liouvillian.jl` | Collected `ρ ↦ AρB` terms, channels, composition, Liouvillian lowering | `Liouvillian`, `terms`, `collapse`, `jump`, `compose`, `harmonics` |
| `engine.jl` | Generic Van Vleck recursion and order scaling | `FloquetExpansion`, `floquet_expansion`, `effective_generator`, `micromotion` |
| `quasienergy.jl` | Symbolic Sambe blocks and harmonic indexing | `QuasienergyOperator`, `harmonic_range` |

The package delegates operator multiplication, adjoints, normal ordering, and coefficient algebra to SecondQuantizedAlgebra. Keep those concerns at that dependency's seam instead of recreating them in this package.

## Representation rules

- A missing Fourier harmonic and a zero harmonic are semantically equivalent.
- A Liouvillian is a collected sum of left/right terms; its sparse dictionary is an implementation detail exposed through `terms`.
- A finite-order effective generator is returned as the raw algebraic truncation. Lindblad form or complete positivity is not inferred or repaired.
- Dissipative quasienergy blocks use the energy-like Floquet-Liouville convention documented in ADR 0008; numerical vectorization remains an adapter concern.

When a proposed change conflicts with this map or an ADR, state the conflict and update the decision record before treating the new behavior as settled.

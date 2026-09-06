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
                                                      ├─ effective_component
                                                      ├─ micromotion
                                                      └─ optional positive_completion
                                                               │
                                                               ▼
                                                        FloquetExpansion
                                                        (completed state)
                                                               │
                                                               ├─ channels
                                                               ├─ kossakowski
                                                               ├─ dissipative_frame
                                                               └─ factorization
```

`PeriodicGenerator{T}` is the shared Fourier boundary. Its component type carries the algebra through addition, commutators, derivatives, antiderivatives, and simplification. The expansion engine should remain independent of whether `T` is `SQA.QAdd` or `Liouvillian`.

Positive completion is a separate opt-in stage after the finite Floquet expansion. It preserves the retained Floquet coefficients and micromotion while replacing the finite effective-generator realization by a selected positive continuation. Completion returns another `FloquetExpansion`; no parallel completed-result wrapper is introduced.

## Module ownership

Listed in `include` order, which is also the dependency order. Every file under `src/` has a row.

| Module | Owns | Public seam |
| --- | --- | --- |
| `FloquetExpansions.jl` | Module wiring, `include` order, the `@reexport` of SecondQuantizedAlgebra, the export list | the exports themselves |
| `periodic_operator.jl` | Fourier harmonics, drive frequency, harmonic calculus, gauges | `PeriodicGenerator`, `harmonics`, `time_average`, `derivative`, `antiderivative` |
| `completion_types.jl` | The completion type lattice, the `PositiveCompletion` field set, provenance types, the completion exceptions | `Completion`, `Uncompleted`, `CompletionAlgorithm`, `Gram`, `Spectral`, `CompletionFactorization` |
| `liouvillian.jl` | Collected `ρ ↦ AρB` terms, physical channel constructors, composition, Liouvillian lowering | `Liouvillian`, `terms`, `collapse`, `jump`, `compose`, `harmonics` |
| `quasienergy.jl` | Symbolic Sambe blocks and harmonic indexing | `QuasienergyOperator`, `harmonic_range` |
| `engine.jl` | Generic Van Vleck recursion, order scaling, `FloquetExpansion`, retained effective/micromotion accessors, the `channels=` keyword | `FloquetExpansion`, `floquet_expansion`, `effective_generator`, `effective_component`, `micromotion` |
| `gksl_coordinates.jl` | Ordered dissipative frames and exact GKSL/Kossakowski coordinate extraction | `DissipativeFrame`, `hamiltonian`, `hamiltonian_component`, `kossakowski`, `kossakowski_component` |
| `completion.jl` | The `positive_completion` dispatch surface and the uncompleted-to-completed state transition | `positive_completion` |
| `gksl_floquet.jl` | GKSL coordinate extraction specialised to a `FloquetExpansion` | `kossakowski`, `kossakowski_component`, `hamiltonian`, `hamiltonian_component` |

**Some representation rules below describe API that does not exist yet on this branch.** `Gram` and `Spectral` are declared in `completion_types.jl` but their algorithms are not implemented here. `positivity_conditions`, `regularity_conditions` and `factorization` are fields of the completion struct, not accessor functions. `channels` is a keyword argument to `floquet_expansion`, not an accessor. There is no `dissipative_frame` function.

So read a rule naming `channels(cp)`, `factorization(cp)`, or the behaviour of `Gram()` and `Spectral()` as the target the completion work is building toward, not as a description of this tree. The rules stay because they constrain that work. Their gate is review, and the code has not reached them yet.

The package delegates operator multiplication, adjoints, normal ordering, and coefficient algebra to SecondQuantizedAlgebra. Keep those concerns at that dependency's seam instead of recreating them in this package.

## Representation rules

- A missing Fourier harmonic and a zero harmonic are semantically equivalent.
- A Liouvillian is a collected sum of left/right terms; its sparse dictionary is an implementation detail exposed through `terms`.
- `Liouvillian` and `PeriodicGenerator` remain algebraic and do not carry dissipative provenance through arbitrary arithmetic.
- The high-level physical `floquet_expansion(...; channels=...)` path may retain internal microscopic channel provenance in the resulting `FloquetExpansion` for later completion.
- A raw finite-order effective generator is returned as the algebraic truncation and is not assumed to be GKSL or completely positive.
- Positive completion is explicit, never implicit in `floquet_expansion`, and does not rewrite retained Floquet coefficients or micromotion.
- Kossakowski coordinates are always relative to an ordered `DissipativeFrame`. Ordering is representation-significant even when two frames span the same subspace.
- `Gram()` is algebraic and must not require Hilbert-space/Liouville-space matrices or symbolic eigendecomposition.
- `Spectral()` is a restricted perturbative spectral/HCM realization and does not define the general completion architecture.
- `channels(cp)` is defined only for completed expansions and must reconstruct `effective_generator(cp)` together with `hamiltonian(cp)`.
- Algorithm-specific intermediate data are exposed only through `factorization(cp)`; the normal completed physical interface is common to all completion algorithms.
- Dissipative quasienergy blocks use the energy-like Floquet-Liouville convention documented in ADR 0008; numerical vectorization remains an adapter concern.

When a proposed change conflicts with this map or an ADR, state the conflict and update the decision record before treating the new behavior as settled.

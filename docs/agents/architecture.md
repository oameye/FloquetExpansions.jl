# Architecture Map

Read this guide when a change crosses module boundaries, changes a symbolic representation, or changes the public expansion/completion model. Read the relevant ADRs for rationale and constraints.

## Data flow

```text
time-dependent Hamiltonian / Liouvillian / PeriodicGenerator
                    │
                    ▼
          PeriodicGenerator{T}
                    │
        ┌───────────┴────────────┐
        ▼                        ▼
 Van Vleck recursion      QuasienergyOperator
        │
        ▼
 FloquetExpansion (raw)
        │
        ├─ effective_generator / effective_component
        ├─ micromotion
        └─ optional positive_completion   [Liouvillian only]
                         │
                         ▼
                FloquetExpansion (completed)
                         │
                         ├─ effective_generator
                         ├─ channels / hamiltonian
                         ├─ kossakowski / dissipative_frame
                         └─ factorization / conditions
```

`PeriodicGenerator{T}` is the shared Fourier boundary. High-level Hamiltonian-plus-channel input is lowered to a Liouvillian periodic generator before recursion; an already constructed `PeriodicGenerator` enters at that boundary directly. Its component type carries the algebra through addition, commutators, derivatives, antiderivatives, and simplification. The expansion engine remains generic over Hamiltonian (`SQA.QAdd`) and Liouvillian components.

Positive completion is an explicit post-processing stage for Liouvillian expansions. It preserves the retained Floquet coefficients and micromotion while replacing only the finite effective-generator realization by a selected positive continuation. Completion returns another `FloquetExpansion`; there is no parallel completed-result wrapper.

## Module ownership

Listed in `src/FloquetExpansions.jl` include order, which is the source dependency order. Every file currently under `src/` has a row. The public-seam column lists names owned by this package; SecondQuantizedAlgebra reexports are wired at the module root.

| Module | Owns | Public seam |
| --- | --- | --- |
| `FloquetExpansions.jl` | Module wiring, include order, SQA reexports/forwards, the export list, and the qualified expert API | the module's exported and `@public` names |
| `periodic_operator.jl` | Fourier harmonics, drive frequency, harmonic calculus, gauges | `Gauge`, `VanVleck`, `PeriodicGenerator`, `harmonics`, `support`, `time_average`, `derivative`, `antiderivative` |
| `completion_types.jl` | Completion state/algorithm types, factorization supertype, completion exceptions, microscopic provenance, completed-state storage | `Completion`, `Uncompleted`, `CompletionAlgorithm`, `Gram`, `Spectral`, `CompletionFactorization`, `CompletionObstruction`, `FractionalJumpOnset` |
| `matrix_series.jl` | Truncated completion scalar/matrix-series algebra, conditions, and graded factor recurrences | internal only |
| `completion_linear_algebra.jl` | Reusable symbolic solve plans, triangular series solves, structured Hermitian congruence elimination, Gram/Feshbach dressing | internal only |
| `liouvillian.jl` | Collected `ρ ↦ AρB` terms, coherent/dissipative constructors, physical channel values, composition, Liouvillian Fourier lowering | `Liouvillian`, `liouvillian`, `terms`, `hamiltonian_action`, `dissipator`, `compose`, `collapse`, `jump`, plus the `harmonics` Liouvillian method |
| `quasienergy.jl` | Symbolic Sambe blocks and harmonic indexing | `QuasienergyOperator`, `harmonic_range` |
| `engine.jl` | Generic Van Vleck recursion, order scaling, `FloquetExpansion`, retained effective/micromotion accessors, high-level physical lowering and microscopic-channel retention | `FloquetExpansion`, `floquet_expansion`, `order`, `effective_generator`, `effective_component`, `micromotion` |
| `gksl_coordinates.jl` | Ordered dissipative frames and exact GKSL/Kossakowski coordinate extraction | `DissipativeFrame`, `hamiltonian`, `hamiltonian_component`, `kossakowski`, `kossakowski_component` |
| `completion_conversion.jl` | Narrow conversion boundary between SQA coefficients and the completion scalar backend | internal only |
| `completion_frame.jl` | Automatic dissipative-frame discovery and independent-direction filtering modulo identity | internal only |
| `gram_completion.jl` | Gram completion data types and single-stratum factor operations | `GramStage`, `GramFactorization` |
| `gram_recursion.jl` | Recursive active/dark onset filtration and the algebraic Gram completion implementation | internal only |
| `spectral_completion.jl` | Restricted perturbative spectral/HCM completion and factorization data | `SpectralFactorization` |
| `completion.jl` | Common completion dispatch, result finalization, owned representation data, retained/coherent caches, and common completed-state accessors | `positive_completion`, `channels`, `dissipative_frame`, `positivity_conditions`, `regularity_conditions`, `factorization` |
| `gksl_floquet.jl` | GKSL/Kossakowski and coherent-Hamiltonian accessors specialized to `FloquetExpansion` | `kossakowski`, `kossakowski_component`, `hamiltonian`, `hamiltonian_component` |

The package delegates operator multiplication, adjoints, normal ordering, and coefficient algebra to SecondQuantizedAlgebra. Keep those concerns at that dependency's seam instead of recreating them here. Completion currently keeps a dedicated symbolic scalar backend behind `completion_conversion.jl`; replacing that backend is a separate architectural change, not a cleanup inside Gram or Spectral code.

## Representation rules

- A missing Fourier harmonic and a zero harmonic are semantically equivalent.
- A Liouvillian is a collected sum of left/right terms; its sparse dictionary is an implementation detail exposed through `terms`.
- `Liouvillian` and `PeriodicGenerator` remain algebraic and do not carry dissipative provenance through arbitrary arithmetic.
- The high-level physical `floquet_expansion(...; channels=...)` path may retain internal microscopic channel provenance in the resulting `FloquetExpansion` for later completion.
- A raw finite-order effective generator is the algebraic truncation and is not assumed to be GKSL or completely positive.
- Positive completion is explicit, never implicit in `floquet_expansion`, is defined only for Liouvillian expansions, and does not rewrite retained Floquet coefficients or micromotion.
- Kossakowski coordinates are relative to an ordered `DissipativeFrame`. Ordering is representation-significant even when two frames span the same subspace.
- Automatic frame discovery is a symbolic convenience frontend whose output arity depends on runtime algebraic independence. The explicit-frame `positive_completion(expansion, algorithm, frame)` methods are the inference-oriented computational core.
- Automatic frame discovery starts from microscopic dissipative directions when provenance is available and appends algebraically generated independent directions deterministically, with independence taken modulo the identity.
- A completed result owns its finalized `DissipativeFrame` and caches the physical retained Kossakowski coefficients and coherent Hamiltonian used to construct its completed generator. Public accessors return defensive copies where mutation could invalidate that owned representation.
- `Gram()` is algebraic and must not require Hilbert-space/Liouville-space matrices, characteristic polynomials, symbolic eigendecomposition, or symbolic matrix square roots.
- `Spectral()` is a restricted perturbative spectral/HCM realization. It requires the leading Kossakowski form to be diagonal in its frame and does not define the general completion architecture.
- `channels(cp)` is defined only for completed expansions and satisfies `liouvillian(hamiltonian(cp); channels=channels(cp)) == effective_generator(cp)`.
- Algorithm-specific intermediate data are exposed through `factorization(cp)`; the normal completed physical interface is common to all completion algorithms.
- Expert API names marked `@public` are stable qualified interfaces, intentionally not widened into the ordinary export list.
- Dissipative quasienergy blocks use the energy-like Floquet-Liouville convention documented in ADR 0008; numerical vectorization remains an adapter concern.

Implementation-performance constraints such as solve-plan reuse and structured Hermitian elimination are owned by [`performance.md`](performance.md); the completion physics and representation contract are owned by ADR 0009.

When a proposed change conflicts with this map or an ADR, state the conflict and update or supersede the decision record before treating the new behavior as settled.

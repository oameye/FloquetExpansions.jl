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
 FloquetExpansion
        │
        ├─ Hamiltonian: canonical HFE
        └─ Liouvillian
             ├─ default: Gram-completed static GKSL generator
             └─ complete_positive=Val(false): raw canonical HFE
                         │
                         ├─ effective_generator / effective_component
                         ├─ micromotion
                         ├─ channels / hamiltonian       [completed only]
                         ├─ kossakowski / dissipative_frame
                         └─ factorization / conditions   [completed only]
```

`PeriodicGenerator{T}` is the shared Fourier boundary. High-level Hamiltonian-plus-channel input is lowered to a Liouvillian periodic generator before recursion; an already constructed `PeriodicGenerator` enters at that boundary directly. Its component type carries the algebra through addition, commutators, derivatives, antiderivatives, and simplification. The expansion engine remains generic over Hamiltonian (`SQA.QAdd`) and Liouvillian components.

For Liouvillian input, the public Hori–Deprit and Bloch/Feshbach Van Vleck selectors apply graded `Gram()` positive completion by default. The raw canonical HFE is selected statically with `complete_positive=Val(false)`. Completion preserves the retained Floquet coefficients and micromotion while replacing only the finite effective-generator realization by a positive continuation. The explicit `positive_completion` API remains available for raw expansions when a fixed frame, `Spectral()`, or an explicit algorithm comparison is required. There is no parallel completed-result wrapper.

## Module ownership

Listed in `src/FloquetExpansions.jl` dependency order. The table emphasizes the core public and policy modules; internal CK reconstruction helpers form the physical open-leg/CPTP branch integrated alongside this generic expansion stack.

| Module | Owns | Public seam |
| --- | --- | --- |
| `FloquetExpansions.jl` | Module wiring, include order, SQA reexports/forwards, the export list, and the qualified expert API | the module's exported and `@public` names |
| `periodic_operator.jl` | Fourier harmonics, drive frequency, harmonic calculus, gauges, algorithm selector types | `Gauge`, `ExpansionAlgorithm`, `HoriDeprit`, `BlochFeshbach`, `VanVleck`, `PeriodicGenerator`, `harmonics`, `support`, `time_average`, `derivative`, `antiderivative` |
| `bloch_projection.jl` | Projection/Feshbach recurrence plans and finite-Fourier projection machinery | internal projection core |
| `bloch_van_vleck.jl` / `bloch_connected_van_vleck.jl` | Bloch/Feshbach canonical Van Vleck reconstruction | internal backend |
| `expansion_algorithms.jl` | Algorithm-specific Hori–Deprit and Bloch/Feshbach Van Vleck dispatch | internal dispatch behind `ExpansionAlgorithm` selectors |
| `cp_algorithm_policy.jl` | Static CP policy wrapper and default Gram-completion dispatch for Liouvillian algorithms | `HoriDeprit`, `BlochFeshbach` constructor policy |
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
- Liouvillian Van Vleck expansion is `Gram()`-completed by default; explicit `complete_positive=Val(false)` selects the raw canonical HFE. Completion does not rewrite retained Floquet coefficients or micromotion.
- Kossakowski coordinates are relative to an ordered `DissipativeFrame`. Ordering is representation-significant even when two frames span the same subspace.
- Raw expansions require an explicit `DissipativeFrame` for GKSL/Kossakowski coordinate extraction. Completed expansions store the finalized frame, so no-frame completed accessors are unambiguous.
- Automatic frame discovery is a symbolic convenience frontend whose output arity depends on runtime algebraic independence. The explicit-frame `positive_completion(expansion, algorithm, frame)` methods are the inference-oriented computational core.
- Automatic frame discovery starts from microscopic dissipative directions when provenance is available and appends algebraically generated independent directions deterministically, with independence taken modulo the identity.
- A completed result owns its finalized `DissipativeFrame` and caches the physical retained Kossakowski coefficients and coherent Hamiltonian used to construct its completed generator. Public accessors return defensive copies where mutation could invalidate that owned representation.
- `Gram()` is algebraic and must not require Hilbert-space/Liouville-space matrices, characteristic polynomials, symbolic eigendecomposition, or symbolic matrix square roots.
- `Spectral()` is a restricted perturbative spectral/HCM realization. The leading Kossakowski form must be diagonal in its frame, and retained corrections must already be diagonal within degenerate leading sectors. Automatic frame discovery does not diagonalize the leading Kossakowski form.
- `channels(cp)` is defined only for completed expansions and satisfies `liouvillian(hamiltonian(cp); channels=channels(cp)) == effective_generator(cp)`.
- Algorithm-specific intermediate data are exposed through `factorization(cp)`; the normal completed physical interface is common to all completion algorithms.
- Native CK/open-leg reconstruction produces a finite CPTP period map through canonical amplitude normalization, output pairing, and physical right normalization. It is distinct from static Gram completion, Floquet gauge choice, and Kraus/channel gauge.
- Expert API names marked `@public` are stable qualified interfaces, intentionally not widened into the ordinary export list.
- Dissipative quasienergy blocks use the energy-like Floquet-Liouville convention documented in ADR 0008; numerical vectorization remains an adapter concern.

Implementation-performance constraints such as solve-plan reuse and structured Hermitian elimination are owned by [`performance.md`](performance.md); the completion physics and representation contract are owned by ADR 0009.

When a proposed change conflicts with this map or an ADR, state the conflict and update or supersede the decision record before treating the new behavior as settled.

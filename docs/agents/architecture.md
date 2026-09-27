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
| `gauges.jl` | Gauges and the expansion-algorithm selectors a gauge carries | `Gauge`, `VanVleck`, `ExpansionAlgorithm`, `HoriDeprit`, `BlochFeshbach` |
| `periodic_operator.jl` | Fourier harmonics, drive frequency, harmonic calculus | `PeriodicGenerator`, `harmonics`, `support`, `time_average`, `derivative`, `antiderivative` |
| `completion_types.jl` | Completion state/algorithm types, factorization supertype, completion exceptions, microscopic provenance, completed-state storage | `Completion`, `Uncompleted`, `CompletionAlgorithm`, `Gram`, `Spectral`, `CompletionFactorization`, `CompletionObstruction`, `FractionalJumpOnset` |
| `matrix_series.jl` | Truncated completion scalar/matrix-series algebra, conditions, and graded factor recurrences | internal only |
| `completion_linear_algebra.jl` | Reusable symbolic solve plans, triangular series solves, structured Hermitian congruence elimination, Gram/Feshbach dressing | internal only |
| `liouvillian.jl` | Collected `ρ ↦ AρB` terms, coherent/dissipative constructors, physical channel values, composition, Liouvillian Fourier lowering | `Liouvillian`, `liouvillian`, `terms`, `hamiltonian_action`, `dissipator`, `compose`, `collapse`, `jump`, plus the `harmonics` Liouvillian method |
| `quasienergy.jl` | Symbolic Sambe blocks and harmonic indexing | `QuasienergyOperator`, `harmonic_range` |
| `floquet_expansion.jl` | `FloquetExpansion`, order scaling, retained effective/micromotion accessors, high-level physical lowering and microscopic-channel retention, the component conventions and input checks shared by the expansion algorithms, and the error for an unimplemented expansion algorithm | `FloquetExpansion`, `floquet_expansion`, `order`, `effective_generator`, `effective_component`, `micromotion` |
| `hori_deprit.jl` | The Hori–Deprit algorithm: the Lie-transform recursion for the van Vleck micromotion and effective generator | internal, selected by `HoriDeprit` |
| `bloch_feshbach/lyndon_words.jl` | Polynomials in the free associative algebra over harmonic letters and their coordinates in the Lyndon commutator basis | internal only |
| `bloch_feshbach/van_vleck_words.jl` | The van Vleck pair over harmonic letters from the support alone: the Bloch recurrence for the wave operator and Bloch effective generator, the connected logarithm and static factor, and the static-factor similarity | internal only |
| `bloch_feshbach/bloch_feshbach.jl` | The Bloch/Feshbach algorithm: compilation of the word-level van Vleck pair to Lyndon commutators, their evaluation on generator components under the Hamiltonian or Liouvillian product and phase conventions, and the entry point | internal, selected by `BlochFeshbach` |
| `gksl_coordinates.jl` | Ordered dissipative frames and exact GKSL/Kossakowski coordinate extraction | `DissipativeFrame`, `hamiltonian`, `hamiltonian_component`, `kossakowski`, `kossakowski_component` |
| `completion_conversion.jl` | Narrow conversion boundary between SQA coefficients and the completion scalar backend | internal only |
| `completion_frame.jl` | Automatic dissipative-frame discovery and independent-direction filtering modulo identity | internal only |
| `gram_completion.jl` | Gram completion data types and single-stratum factor operations | `GramStage`, `GramFactorization` |
| `gram_recursion.jl` | Recursive active/dark onset filtration and the algebraic Gram completion implementation | internal only |
| `spectral_completion.jl` | Restricted perturbative spectral/HCM completion and factorization data | `SpectralFactorization` |
| `completion.jl` | Common completion dispatch, result finalization, owned representation data, retained/coherent caches, and common completed-state accessors | `positive_completion`, `channels`, `dissipative_frame`, `positivity_conditions`, `regularity_conditions`, `factorization` |
| `gksl_floquet.jl` | GKSL/Kossakowski and coherent-Hamiltonian accessors specialized to `FloquetExpansion` | `kossakowski`, `kossakowski_component`, `hamiltonian`, `hamiltonian_component` |

`@reexport` forwards only the names SecondQuantizedAlgebra exports, so the module root imports and exports `expim`, `exponential_form`, and `trigonometric_form`, which SecondQuantizedAlgebra marks `@public` without exporting.

The package delegates operator multiplication, adjoints, normal ordering, and coefficient algebra to SecondQuantizedAlgebra. Keep those concerns at that dependency's seam instead of recreating them here. Completion currently keeps a dedicated symbolic scalar backend behind `completion_conversion.jl`; replacing that backend is a separate architectural change, not a cleanup inside Gram or Spectral code.

## Representation rules

- A missing Fourier harmonic and a zero harmonic are semantically equivalent.
- `iszero` on a Symbolics `BasicSymbolic` builds the symbolic equation `0 == 0` rather than returning a `Bool`, so Fourier lowering tests symbolic zeros structurally.
- `expim(arg)` is ``e^{+i\,\mathrm{arg}}`` while the package Fourier convention is ``e^{-im\omega t}``, so Fourier lowering reads the harmonic label as the negated coefficient of ``\omega t`` and rejects a phase that is not linear in ``\omega t``.
- A Liouvillian is a collected sum of left/right terms; its sparse dictionary is an implementation detail exposed through `terms`.
- `Liouvillian` and `PeriodicGenerator` remain algebraic and do not carry dissipative provenance through arbitrary arithmetic. The provenance types live in `completion_types.jl` for that reason.
- The LaTeX display of a Liouvillian embeds SecondQuantizedAlgebra and Symbolics `text/latex` output inside a larger expression. Both emit self-delimited math, and a renderer strips only the outermost delimiter pair, so the display strips their delimiters through the `show` path; Latexify's environment entry point has no method for some symbolic types.
- The high-level physical `floquet_expansion(...; channels=...)` path may retain internal microscopic channel provenance in the resulting `FloquetExpansion` for later completion.
- A raw finite-order effective generator is the algebraic truncation and is not assumed to be GKSL or completely positive.
- Bloch/Feshbach computes the van Vleck pair from the harmonic support and order alone, rebuilt on every call, in the free associative algebra over harmonic letters, with exact `Rational{Int}` coefficients, and evaluates only Lyndon commutators of the generator components. This relies on the connected logarithm and the van Vleck effective generator both being Lie series in the harmonics; the Lyndon decomposition throws when a word series is not one.
- A `ComponentConvention` fixes, per component type, the component product and the unit phase ``c`` with ``\mathcal{G} = cX``: ``c = -i`` for a Hamiltonian, ``c = 1`` for a Liouvillian. Hori–Deprit dresses with the phase ``-c``. Bloch/Feshbach evaluates words with ``1/m`` in place of ``i/m``, so it restores the phase ``(ic)^n`` on order-``n`` coefficients and a further ``\bar c = c^{-1}`` on the micromotion.
- Positive completion is explicit, never implicit in `floquet_expansion`, is defined only for Liouvillian expansions, and does not rewrite retained Floquet coefficients or micromotion.
- Kossakowski coordinates are relative to an ordered `DissipativeFrame`. Ordering is representation-significant even when two frames span the same subspace.
- Raw expansions require an explicit `DissipativeFrame` for GKSL/Kossakowski coordinate extraction. Completed expansions store the finalized frame, so no-frame completed accessors are unambiguous.
- Automatic frame discovery is a symbolic convenience frontend whose output arity depends on runtime algebraic independence. The explicit-frame `positive_completion(expansion, algorithm, frame)` methods are the inference-oriented computational core.
- Automatic frame discovery starts from microscopic dissipative directions when provenance is available and appends algebraically generated independent directions deterministically, with independence taken modulo the identity.
- The public `DissipativeFrame` constructor rejects an empty frame, but automatic frame discovery returns one for an expansion without dissipative terms, so a coherent Liouvillian completes with no channels and a `0 × 0` Kossakowski matrix. GKSL extraction in the empty frame still rejects any two-sided residual.
- A completed result owns its finalized `DissipativeFrame` and caches the physical retained Kossakowski coefficients and coherent Hamiltonian used to construct its completed generator. Completion copies the caller's frame, so mutating that frame afterwards cannot invalidate the cached coordinates. Public accessors return defensive copies where mutation could invalidate that owned representation.
- Retained GKSL data in a fixed frame hold the coherent part already reattached to physical inverse-drive powers, while the dissipative coefficients stay order-resolved until completion finalization converts them to their physical series.
- `Gram()` is algebraic and must not require Hilbert-space/Liouville-space matrices, characteristic polynomials, symbolic eigendecomposition, or symbolic matrix square roots.
- `Spectral()` is a restricted perturbative spectral/HCM realization. The leading Kossakowski form must be diagonal in its frame, and retained corrections must already be diagonal within degenerate leading sectors. Automatic frame discovery does not diagonalize the leading Kossakowski form.
- `channels(cp)` is defined only for completed expansions and satisfies `liouvillian(hamiltonian(cp); channels=channels(cp)) == effective_generator(cp)`.
- Algorithm-specific intermediate data are exposed through `factorization(cp)`; the normal completed physical interface is common to all completion algorithms.
- Expert API names marked `@public` are stable qualified interfaces, intentionally not widened into the ordinary export list.
- Dissipative quasienergy blocks use the energy-like Floquet-Liouville convention documented in ADR 0008; numerical vectorization remains an adapter concern.

Implementation-performance constraints such as solve-plan reuse and structured Hermitian elimination are owned by [`performance.md`](performance.md); the completion physics and representation contract are owned by ADR 0009.

When a proposed change conflicts with this map or an ADR, state the conflict and update or supersede the decision record before treating the new behavior as settled.

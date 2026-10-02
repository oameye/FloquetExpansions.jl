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

`src/` is one module split into four layer folders plus a shared word algebra. `src/FloquetExpansions.jl` includes them in the order below, which is the source dependency order.

| Folder | Layer |
| --- | --- |
| `generators/` | Periodic generators: Fourier harmonics and harmonic calculus, the `Liouvillian` map algebra, dissipative channels, and Sambe blocks. Knows nothing of gauges, expansions, or completion. |
| `words/` | The free associative algebra over harmonic letters and the word-level van Vleck pair. Generic over its coefficient type and independent of the generator component. |
| `expansion/` | Gauges, the `FloquetExpansion` result, and the expansion algorithms that fill it. |
| `gksl/` | Ordered dissipative frames and exact GKSL/Kossakowski coordinates, on a `Liouvillian` and on any `FloquetExpansion`, plus the shared floating-point Gram geometry. |
| `native/` | The GKSL normal form: state solved inside the BF/HD recurrence rather than by completion, its exact cutoff-free algebra, and the `GKSLNormalForm` entry point. Does not depend on `completion/`. |
| `completion/` | Positive completion: algorithm selectors, the completed state, its accessors, and the Gram and Spectral realizations. `completion/backend/` holds the dedicated completion scalar backend. |

**Layering rule.** A file names a type from its own folder or an earlier one only. A call into a later folder resolves at run time and is allowed, but a signature, a field type, or a type parameter never refers forward. A change that needs a forward type reference moves the type to the earliest folder that uses it, as `Completion` sits in `expansion/` and the provenance types in `generators/`.

Every file currently under `src/` has a row. The public-seam column lists names owned by this package; SecondQuantizedAlgebra reexports are wired at the module root.

| File | Owns | Public seam |
| --- | --- | --- |
| `FloquetExpansions.jl` | Module wiring, include order, SQA reexports/forwards, the export list, and the qualified expert API | the module's exported and `@public` names |
| `generators/periodic_generator.jl` | Fourier harmonics, drive frequency, gauge-free harmonic calculus | `PeriodicGenerator`, `harmonics`, `support`, `time_average`, `derivative` |
| `generators/liouvillian.jl` | Collected `ρ ↦ AρB` terms, coherent and dissipator constructors, composition, Liouvillian Fourier lowering | `Liouvillian`, `terms`, `hamiltonian_action`, `dissipator`, `compose`, plus the `harmonics` Liouvillian method |
| `generators/channels.jl` | Physical channel values and their displays, jump-rate validation, microscopic provenance types and their construction, the Hamiltonian-plus-channels constructor | `liouvillian`, `collapse`, `jump` |
| `generators/quasienergy.jl` | Symbolic Sambe blocks and harmonic indexing | `QuasienergyOperator`, `harmonic_range` |
| `words/lyndon_words.jl` | Polynomials in the free associative algebra over harmonic letters and their coordinates in the Lyndon commutator basis | internal only |
| `words/van_vleck_words.jl` | The van Vleck pair over harmonic letters from the support alone: the Bloch recurrence for the wave operator and Bloch effective generator, the connected logarithm and static factor, and the static-factor similarity | internal only |
| `expansion/gauges.jl` | Gauges, the expansion-algorithm selectors a gauge carries, and the gauge-fixed antiderivative | `Gauge`, `VanVleck`, `ExpansionAlgorithm`, `HoriDeprit`, `BlochFeshbach`, `antiderivative` |
| `expansion/floquet_expansion.jl` | The completion-state supertype and `Uncompleted`, `FloquetExpansion`, order scaling, retained effective/micromotion accessors, high-level physical lowering and microscopic-channel retention, the component conventions and input checks shared by the expansion algorithms, and the error for an unimplemented expansion algorithm | `FloquetExpansion`, `floquet_expansion`, `order`, `effective_generator`, `effective_component`, `micromotion`, `Completion`, `Uncompleted` |
| `expansion/hori_deprit.jl` | The Hori–Deprit algorithm: the Lie-transform recursion for the van Vleck micromotion and effective generator | internal, selected by `HoriDeprit` |
| `expansion/bloch_feshbach.jl` | The Bloch/Feshbach algorithm: compilation of the word-level van Vleck pair to Lyndon commutators, their evaluation on generator components under the Hamiltonian or Liouvillian product and phase conventions, and the entry point | internal, selected by `BlochFeshbach` |
| `gksl/dissipative_frame.jl` | Ordered dissipative frames, the GKSL coordinate error, and the exact coefficient linear algebra that builds and inverts frame coordinates | `DissipativeFrame` |
| `gksl/coordinates.jl` | Exact GKSL/Kossakowski coordinate extraction from a `Liouvillian` in a frame | `hamiltonian`, `kossakowski` on a `Liouvillian` |
| `gksl/floquet.jl` | GKSL/Kossakowski and coherent-Hamiltonian accessors on a `FloquetExpansion` in an explicit frame, for any completion state | `kossakowski`, `kossakowski_component`, `hamiltonian`, `hamiltonian_component` |
| `gksl/gram_kernel.jl` | Floating-point Gram geometry of fixed amplitude columns and Hermitian affine slices: active/dark frames, tangent-Gram lifts, PSD factors, and the affine PSD section | internal only |
| `native/graded_channels.jl` | Graded jump amplitudes with half-integer onset, their active flag, known products ``K_n^<``, Gram coefficients, and prefix-complete storage of corrections and births | internal only |
| `native/static_solve.jl` | One native static slot in Kossakowski coordinates: dark quotient, affine PSD gauge section, newborn channels, bright tangent lift, and reconstruction check. The superoperator representation stays with the caller | internal only |
| `native/exact_static_solve.jl` | The exact counterpart of `native_static_solve`, generic over the scalar field: facial reduction of the dark positivity problem, the least-squares section with a minimum-norm tie-break and a row-cancellation fallback, rate-weighted ``LDL^\dagger`` births, and the closed-form no-bright-rotation lift. `NativePositivityError` reports a gauge family without a PSD lift | internal only |
| `native/static_step.jl` | The generic native static step over a `NativeRepresentation` (Kossakowski map, Hamiltonian read-off, GKSL assembly, gauge directions) and the `HomologicalInverse` seam: an inverse returns a static gauge whose commutator with ``\bar{\mathcal L}`` cancels the dark target exactly on the sectors where ``\operatorname{ad}_{\bar{\mathcal L}}`` is regular, and the gauge directions spanning the remaining sector, which keeps the affine PSD section. `NoHomologicalInverse` has no regular sector | internal only |
| `native/dense_representation.jl` | The finite-dimensional dense Liouville representation: Hilbert–Schmidt traceless basis, Kossakowski and Hamiltonian extraction, GKSL assembly, and the Hermiticity- and trace-preserving static gauge algebra, built in the constructor | internal only |
| `native/exact_representation.jl` | The exact dense Liouville representation with an explicit diagonal Hilbert metric, so rescaled Fock and spin bases stay rational | internal only |
| `native/charge_grading.jl` | Detection of the U(1) charge lattice of a dense averaged generator and its leading channels, and the resolvent section on regular charged blocks | internal only |
| `native/operator_algebra.jl` | The exact cutoff-free operator algebra: boson, phase-space, and finite-level sites, normal-ordered monomials, operators, and left/right superoperators | internal only |
| `native/algebraic_representation.jl` | The native representation over that algebra: a finite monomial frame, exact Kossakowski and Hamiltonian read-off, polynomial static gauge families by degree, and a polynomial charged inverse. `NativeFrameError` reports a frame that is too small | internal only |
| `native/sqa_lowering.jl` | Lowering of SQA operators and periodic Liouvillians into the exact algebra with exact coefficients, and lifting of results back to SQA expressions | internal only |
| `native/recurrence.jl` | The arbitrary-order native recurrences: Bloch–Feshbach and Hori–Deprit over harmonic series of representation elements, with the intrinsic total-kick offset, one native static step per order, and graded-channel bookkeeping. Selected by the existing `BlochFeshbach` and `HoriDeprit` expansion algorithms | internal only |
| `native/driver.jl` | The exact native driver: leading channels from the averaged Kossakowski factor, frame growth, and Int128-to-BigInt promotion | internal only |
| `native/realization.jl` | `floquet_expansion` for `GKSLNormalForm`: lifts effective components, micromotion, Hamiltonian, and channels into a `FloquetExpansion` carrying a `NativeRealization` | `GKSLNormalForm` via `floquet_expansion` |
| `completion/types.jl` | Completion algorithm selectors, factorization supertype, completion exceptions, retained GKSL data, completed-state storage | `CompletionAlgorithm`, `Gram`, `Spectral`, `CompletionFactorization`, `CompletionObstruction`, `FractionalJumpOnset` |
| `completion/backend/matrix_series.jl` | Truncated completion scalar/matrix-series algebra, conditions, and graded factor recurrences | internal only |
| `completion/backend/linear_algebra.jl` | Reusable symbolic solve plans, triangular series solves, structured Hermitian congruence elimination, Gram/Feshbach dressing | internal only |
| `completion/backend/conversion.jl` | Narrow conversion boundary between SQA coefficients and the completion scalar backend | internal only |
| `completion/frame_discovery.jl` | Automatic dissipative-frame discovery and independent-direction filtering modulo identity | internal only |
| `completion/accessors.jl` | Common completed-state accessors, shared by every completion algorithm, and the completed effective generator | `channels`, `dissipative_frame`, `positivity_conditions`, `regularity_conditions`, `factorization`, plus the no-frame `kossakowski`, `kossakowski_component`, `hamiltonian`, and `effective_generator` methods for a completed expansion |
| `completion/gram/factorization.jl` | Gram completion data types and single-stratum factor operations | `GramStage`, `GramFactorization` |
| `completion/gram/recursion.jl` | Recursive active/dark onset filtration and the algebraic Gram completion implementation | internal only |
| `completion/spectral.jl` | Restricted perturbative spectral/HCM completion and factorization data | `SpectralFactorization` |
| `completion/positive_completion.jl` | Common completion dispatch, result finalization, owned representation data, and retained/coherent caches | `positive_completion` |

### Where new work goes

- A new expansion algorithm for an existing gauge is one file in `expansion/` plus its selector in `expansion/gauges.jl`. A new gauge adds its type and its `antiderivative` method to `expansion/gauges.jl`.
- Algebra over harmonic words that another algorithm could reuse belongs in `words/`, not in the algorithm file that first needs it.
- An open-system construction that produces a GKSL generator directly, without passing through `positive_completion`, gets its own folder after `gksl/` and does not depend on `completion/`.
- Replacing the completion scalar backend replaces `completion/backend/` and nothing outside it.
- Multi-frequency harmonics change `generators/periodic_generator.jl` and the harmonic letters in `words/`; the layers above see harmonics only through `PeriodicGenerator`.

`@reexport` forwards only the names SecondQuantizedAlgebra exports, so the module root imports and exports `expim`, `exponential_form`, and `trigonometric_form`, which SecondQuantizedAlgebra marks `@public` without exporting.

The package delegates operator multiplication, adjoints, normal ordering, and coefficient algebra to SecondQuantizedAlgebra. Keep those concerns at that dependency's seam instead of recreating them here. Completion currently keeps a dedicated symbolic scalar backend behind `completion/backend/conversion.jl`; replacing that backend is a separate architectural change confined to `completion/backend/`, not a cleanup inside Gram or Spectral code.

## Representation rules

- A missing Fourier harmonic and a zero harmonic are semantically equivalent.
- `iszero` on a Symbolics `BasicSymbolic` builds the symbolic equation `0 == 0` rather than returning a `Bool`, so Fourier lowering tests symbolic zeros structurally.
- `expim(arg)` is ``e^{+i\,\mathrm{arg}}`` while the package Fourier convention is ``e^{-im\omega t}``, so Fourier lowering reads the harmonic label as the negated coefficient of ``\omega t`` and rejects a phase that is not linear in ``\omega t``.
- A Liouvillian is a collected sum of left/right terms keyed by pairs of unit monomials (ADR 0007); its sparse dictionary is an implementation detail exposed through `terms`.
- `Liouvillian` and `PeriodicGenerator` remain algebraic and do not carry dissipative provenance through arbitrary arithmetic. The provenance types live in `generators/channels.jl`, beside the channels they record, for that reason.
- The LaTeX display of a Liouvillian embeds SecondQuantizedAlgebra and Symbolics `text/latex` output inside a larger expression. Both emit self-delimited math, and a renderer strips only the outermost delimiter pair, so the display strips their delimiters through the `show` path; Latexify's environment entry point has no method for some symbolic types.
- The high-level physical `floquet_expansion(...; channels=...)` path may retain internal microscopic channel provenance in the resulting `FloquetExpansion` for later completion.
- A raw finite-order effective generator is the algebraic truncation and is not assumed to be GKSL or completely positive.
- Bloch/Feshbach computes the van Vleck pair from the harmonic support and order alone, rebuilt on every call, in the free associative algebra over harmonic letters, with exact `Rational{Int}` coefficients, and evaluates only Lyndon commutators of the generator components. This relies on the connected logarithm and the van Vleck effective generator both being Lie series in the harmonics; the Lyndon decomposition throws when a word series is not one.
- A `ComponentConvention` fixes, per component type, the component product and the unit phase ``c`` with ``\mathcal{G} = cX``: ``c = -i`` for a Hamiltonian, ``c = 1`` for a Liouvillian. Hori–Deprit dresses with the phase ``-c``. Bloch/Feshbach evaluates words with ``1/m`` in place of ``i/m``, so it restores the phase ``(ic)^n`` on order-``n`` coefficients and a further ``\bar c = c^{-1}`` on the micromotion.
- Positive completion is explicit, never implicit in `floquet_expansion`, is defined only for Liouvillian expansions, and does not rewrite retained Floquet coefficients or micromotion.
- Kossakowski coordinates are relative to an ordered `DissipativeFrame`. Ordering is representation-significant even when two frames span the same subspace.
- Raw expansions require an explicit `DissipativeFrame` for GKSL/Kossakowski coordinate extraction. Completed expansions store the finalized frame, so no-frame completed accessors are unambiguous.
- Automatic frame discovery is a symbolic convenience frontend whose output arity depends on runtime algebraic independence. The explicit-frame `positive_completion(expansion, algorithm, frame)` methods are the inference-oriented computational core.
- Automatic frame discovery starts from microscopic dissipative directions when provenance is available and appends algebraically generated independent directions deterministically, with independence taken modulo the identity. The microscopic directions are the Fourier harmonics of each channel operator, in user channel order and with the time average of each channel first, so a time-dependent channel never enters the frame with its phases.
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

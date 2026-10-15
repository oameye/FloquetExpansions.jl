# Architecture Map

Read this guide when a change crosses module boundaries, changes a symbolic representation, or changes the public expansion/completion model. Read the relevant ADRs for rationale and constraints.

## Data flow

```text
time-dependent Hamiltonian / Lindbladian / SuperOperator / PeriodicGenerator
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
                         ├─ Lindbladian / channels / hamiltonian
                         ├─ kossakowski / jump_basis
                         └─ factorization / conditions
```

`PeriodicGenerator{T}` is the shared Fourier boundary. A SecondQuantizedAlgebra `Lindbladian` or `SuperOperator` is lowered to a `PeriodicGenerator{SQA.SuperOperator}` before recursion; an already constructed `PeriodicGenerator` enters at that boundary directly. Its component type carries the algebra through addition, commutators, derivatives, antiderivatives, and simplification. The expansion engine remains generic over Hamiltonian (`SQA.QAdd`) and Liouvillian (`SQA.SuperOperator`) components. FloquetExpansions defines no map type of its own (ADR 0014).

Positive completion is an explicit post-processing stage for Liouvillian expansions. It preserves the retained Floquet coefficients and micromotion while replacing only the finite effective-generator realization by a selected positive continuation. Completion returns another `FloquetExpansion`; there is no parallel completed-result wrapper.

## Module ownership

`src/` is one module split into seven layer folders. `src/FloquetExpansions.jl` includes them in the order below, which is the source dependency order.

| Folder | Layer |
| --- | --- |
| `generators/` | Periodic generators: Fourier harmonics and harmonic calculus, Fourier lowering of SQA maps, microscopic channel provenance, and Sambe blocks. Knows nothing of gauges, expansions, or completion. |
| `words/` | The free associative algebra over harmonic letters and the word-level van Vleck pair. Generic over its coefficient type and independent of the generator component. |
| `expansion/` | Gauges, the `FloquetExpansion` result, and the expansion algorithms that fill it. |
| `gksl/` | GKSL/Kossakowski accessors on a `FloquetExpansion` in an SQA `JumpBasis`, and the coefficient-matrix helpers completion needs. |
| `native/` | The GKSL normal form: state solved inside the BF/HD recurrence rather than by completion, its exact cutoff-free algebra, and the `GKSLNormalForm` entry point. Does not depend on `completion/`. |
| `kraus/` | The Floquet–Kraus period map: the no-jump Floquet frame, exact jump-time kernels on the simplex and their Gram matrix, the period-map superoperator, the `period_map` constructor and its `floquet_expansion` method. Reuses the `native/` algebra and recurrence; does not depend on `completion/`. |
| `completion/` | Positive completion: algorithm selectors, the completed state, its accessors, and the Gram and Spectral realizations. `completion/backend/` holds the dedicated completion scalar backend. |

**Layering rule.** A file names a type from its own folder or an earlier one only. A call into a later folder resolves at run time and is allowed, but a signature, a field type, or a type parameter never refers forward. A change that needs a forward type reference moves the type to the earliest folder that uses it, as `Completion` sits in `expansion/` and the provenance types in `generators/`.

Every file currently under `src/` has a row. The public-seam column lists names owned by this package; SecondQuantizedAlgebra reexports are wired at the module root.

| File | Owns | Public seam |
| --- | --- | --- |
| `FloquetExpansions.jl` | Module wiring, include order, SQA reexports/forwards, the export list, and the qualified expert API | the module's exported and `@public` names |
| `generators/periodic_generator.jl` | Fourier harmonics, drive frequency, gauge-free harmonic calculus | `PeriodicGenerator`, `harmonics`, `support`, `time_average`, `derivative` |
| `generators/superoperator.jl` | Sandwich iteration over an SQA `SuperOperator`, completeness canonicalization of produced maps, the canonical composition and commutator used by the expansion engine, Fourier lowering of maps | the `harmonics` methods for `Lindbladian` and `SuperOperator` |
| `generators/channels.jl` | Microscopic provenance types and their construction from the channels of a `Lindbladian`, channel-rate validation, assembly of a `Lindbladian` from a Hamiltonian and `jump => rate` pairs | internal only |
| `generators/quasienergy.jl` | Symbolic Sambe blocks and harmonic indexing | `QuasienergyOperator`, `harmonic_range` |
| `words/lyndon_words.jl` | Polynomials in the free associative algebra over harmonic letters and their coordinates in the Lyndon commutator basis | internal only |
| `words/van_vleck_words.jl` | The van Vleck pair over harmonic letters from the support alone: the Bloch recurrence for the wave operator and Bloch effective generator, the connected logarithm and static factor, and the static-factor similarity | internal only |
| `expansion/gauges.jl` | Gauges, the expansion-algorithm selectors a gauge carries, and the gauge-fixed antiderivative | `Gauge`, `VanVleck`, `GKSLNormalForm`, `ExpansionAlgorithm`, `HoriDeprit`, `BlochFeshbach`, `antiderivative` |
| `expansion/floquet_expansion.jl` | The completion-state supertype and `Uncompleted`, `FloquetExpansion`, order scaling, retained effective/micromotion accessors, high-level physical lowering and microscopic-channel retention, the component conventions and input checks shared by the expansion algorithms, and the error for an unimplemented expansion algorithm | `FloquetExpansion`, `floquet_expansion`, `order`, `effective_generator`, `effective_component`, `micromotion`, `Completion`, `Uncompleted` |
| `expansion/hori_deprit.jl` | The Hori–Deprit algorithm: the Lie-transform recursion for the van Vleck micromotion and effective generator | internal, selected by `HoriDeprit` |
| `expansion/bloch_feshbach.jl` | The Bloch/Feshbach algorithm: compilation of the word-level van Vleck pair to Lyndon commutators, their evaluation on generator components under the Hamiltonian or Liouvillian product and phase conventions, and the entry point | internal, selected by `BlochFeshbach` |
| `gksl/coefficients.jl` | Exact coefficient-matrix helpers shared by completion: products, adjoints, simplification, and the Hermiticity check of a completed Kossakowski matrix | internal only |
| `gksl/floquet.jl` | GKSL/Kossakowski and coherent-Hamiltonian accessors on a `FloquetExpansion` in an explicit jump basis, for any completion state | `kossakowski`, `kossakowski_component`, `hamiltonian`, `hamiltonian_component` |
| `native/graded_channels.jl` | Graded jump amplitudes with half-integer onset, their active flag, known products ``K_n^<``, Gram coefficients, and prefix-complete storage of corrections and births | internal only |
| `native/exact_static_solve.jl` | One native static slot in Kossakowski coordinates, generic over the exact scalar field: dark quotient, facial reduction of the dark positivity problem, the least-squares section with a minimum-norm tie-break and a row-cancellation fallback, rate-weighted ``LDL^\dagger`` births, and the bright tangent lift. `NativePositivityError` reports a gauge family without a PSD lift | internal only |
| `native/static_step.jl` | The generic native static step over a `NativeRepresentation` (Kossakowski map, Hamiltonian read-off, GKSL assembly, gauge families) | internal only |
| `native/operator_algebra.jl` | The exact cutoff-free operator algebra: boson, phase-space, and finite-level sites, normal-ordered monomials, operators, and left/right superoperators | internal only |
| `native/algebraic_representation.jl` | The native representation over that algebra: a finite monomial frame, exact Kossakowski and Hamiltonian read-off, and polynomial static gauge families by degree. `NativeFrameError` reports a frame that is too small | internal only |
| `native/sqa_lowering.jl` | Lowering of SQA operators and periodic Liouvillians into the exact algebra with exact coefficients, and lifting of results back to SQA expressions | internal only |
| `native/recurrence.jl` | The arbitrary-order Bloch–Feshbach and Hori–Deprit recurrences over harmonic series of any associative algebra element, with the intrinsic total-kick offset and a pluggable static slot: the GKSL slot runs one native static step per order with graded-channel bookkeeping, the van Vleck slot sets the static gauge to zero. Selected by the existing `BlochFeshbach` and `HoriDeprit` expansion algorithms | internal only |
| `native/driver.jl` | The exact native driver: leading channels from the averaged Kossakowski factor, frame growth, Int128-to-BigInt promotion, and the coefficient lift that maps exact scalars to SQA coefficients | internal only |
| `native/parameter_field.jl` | Exact symbolic parameters: `ParameterScalar` over Nemo's rational-function field, sign decisions fixed by the declared sign region, the parameter scan and the lift of rational functions back to SQA coefficients | `negative` of `GKSLNormalForm` via `floquet_expansion` |
| `native/realization.jl` | `floquet_expansion` for `GKSLNormalForm`: lifts effective components, micromotion, Hamiltonian, and channels into a `FloquetExpansion` carrying a `NativeRealization` | `GKSLNormalForm` via `floquet_expansion` |
| `kraus/time_kernel.jl` | Operator-valued kernels on the jump-time simplex, graded by inverse-frequency order: static drift between jumps, transported jump insertion, and endpoint conjugation. Exact simplex integrals of the kernel monomials and their Gram matrix | internal only |
| `kraus/no_jump_frame.jl` | The van Vleck frame of the no-jump generator through the shared native recurrence with a zero static slot: effective generator, kick and inverse kick, transported jump operators, and powers of the effective generator | internal only |
| `kraus/period_map.jl` | `PeriodMap` and `PeriodMapExpansion`, the lowering of a `Lindbladian` and its channels into the exact algebra over the rational functions of ``\pi`` and the parameters, the Gram-weighted superoperator sectors and the continuous Kraus amplitudes by jump number, the exact trace-preservation certificate, and the lift back to SQA expressions | `PeriodMap`, `PeriodMapExpansion`, `period_map`, the `floquet_expansion` method on a `PeriodMap`, `superoperator`, `trace_defect`, `JumpAmplitude`, `jump_amplitudes`, `weight`, `amplitude`, `jump_times`, the `order` method on a `PeriodMapExpansion` |
| `completion/types.jl` | Completion algorithm selectors, factorization supertype, completion exceptions, retained GKSL data, completed-state storage | `CompletionAlgorithm`, `Gram`, `Spectral`, `CompletionFactorization`, `CompletionObstruction`, `FractionalJumpOnset` |
| `completion/backend/matrix_series.jl` | Truncated completion scalar/matrix-series algebra, conditions, and graded factor recurrences | internal only |
| `completion/backend/linear_algebra.jl` | Reusable symbolic solve plans, triangular series solves, structured Hermitian congruence elimination, Gram/Feshbach dressing | internal only |
| `completion/backend/conversion.jl` | Narrow conversion boundary between SQA coefficients and the completion scalar backend | internal only |
| `completion/basis_discovery.jl` | Automatic jump-basis discovery from provenance seeds and `support_basis`, filtered by SQA's `independent_basis` | internal only |
| `completion/accessors.jl` | Common completed-state accessors, shared by every completion algorithm, and the completed effective generator | the `Lindbladian` constructor method, `channels`, `jump_basis`, `positivity_conditions`, `regularity_conditions`, `factorization`, plus the no-basis `kossakowski`, `kossakowski_component`, `hamiltonian`, and `effective_generator` methods for a completed expansion |
| `completion/gram/factorization.jl` | Gram completion data types and single-stratum factor operations | `GramStage`, `GramFactorization` |
| `completion/gram/recursion.jl` | Recursive active/dark onset filtration and the algebraic Gram completion implementation | internal only |
| `completion/spectral.jl` | Restricted perturbative spectral/HCM completion and factorization data | `SpectralFactorization` |
| `completion/positive_completion.jl` | Common completion dispatch, result finalization, owned representation data, and retained/coherent caches | `positive_completion` |

### Where new work goes

- A new expansion algorithm for an existing gauge is one file in `expansion/` plus its selector in `expansion/gauges.jl`. A new gauge adds its type and its `antiderivative` method to `expansion/gauges.jl`.
- Algebra over harmonic words that another algorithm could reuse belongs in `words/`, not in the algorithm file that first needs it.
- An open-system construction that produces a GKSL generator directly, without passing through `positive_completion`, belongs in `native/` and does not depend on `completion/`.
- A map-level construction that returns a completely positive map rather than a generator belongs in `kraus/`.
- Replacing the completion scalar backend replaces `completion/backend/` and nothing outside it.
- Multi-frequency harmonics change `generators/periodic_generator.jl` and the harmonic letters in `words/`; the layers above see harmonics only through `PeriodicGenerator`.

`@reexport` forwards only the names SecondQuantizedAlgebra exports, so the module root imports and exports `expim`, `exponential_form`, and `trigonometric_form`, which SecondQuantizedAlgebra marks `@public` without exporting.

The package delegates operator multiplication, adjoints, normal ordering, and coefficient algebra to SecondQuantizedAlgebra. Keep those concerns at that dependency's seam instead of recreating them here. Completion currently keeps a dedicated symbolic scalar backend behind `completion/backend/conversion.jl`; replacing that backend is a separate architectural change confined to `completion/backend/`, not a cleanup inside Gram or Spectral code.

## Representation rules

- A missing Fourier harmonic and a zero harmonic are semantically equivalent.
- `iszero` on a Symbolics `BasicSymbolic` builds the symbolic equation `0 == 0` rather than returning a `Bool`, so Fourier lowering tests symbolic zeros structurally.
- `expim(arg)` is ``e^{+i\,\mathrm{arg}}`` while the package Fourier convention is ``e^{-im\omega t}``, so Fourier lowering reads the harmonic label as the negated coefficient of ``\omega t`` and rejects a phase that is not linear in ``\omega t``.
- A Liouvillian component is an SQA `SuperOperator`. Every map FloquetExpansions produces is passed through `expand_completeness`, so `iszero` and `==` are map equality (ADR 0014). Sandwiches are read through `left_operator` and `right_operator`, and maps are assembled with `SuperBuilder`.
- `SuperOperator` and `PeriodicGenerator` remain algebraic and do not carry dissipative provenance through arbitrary arithmetic. The provenance types live in `generators/channels.jl` for that reason.
- `floquet_expansion` and `period_map` on a `Lindbladian` retain internal microscopic channel provenance in the result for later completion. Only diagonal, unsummed channels with provably real rates are accepted; anything else is passed as `SuperOperator(L)`.
- A raw finite-order effective generator is the algebraic truncation and is not assumed to be GKSL or completely positive.
- Bloch/Feshbach computes the van Vleck pair from the harmonic support and order alone, rebuilt on every call, in the free associative algebra over harmonic letters, with exact `Rational{Int}` coefficients, and evaluates only Lyndon commutators of the generator components. This relies on the connected logarithm and the van Vleck effective generator both being Lie series in the harmonics; the Lyndon decomposition throws when a word series is not one.
- A `ComponentConvention` fixes, per component type, the component product and the unit phase ``c`` with ``\mathcal{G} = cX``: ``c = -i`` for a Hamiltonian, ``c = 1`` for a Liouvillian. Hori–Deprit dresses with the phase ``-c``. Bloch/Feshbach evaluates words with ``1/m`` in place of ``i/m``, so it restores the phase ``(ic)^n`` on order-``n`` coefficients and a further ``\bar c = c^{-1}`` on the micromotion.
- The native layer accepts exact `Integer` and `Rational` coefficients and polynomials in at most eight symbolic parameters, which it expands over the rational-function field (ADR 0012). An order with no static gauge of polynomial degree at most three raises an `ArgumentError`. A `GKSLNormalForm` expansion is already of GKSL form, so `positive_completion` refuses it.
- A period map is exact over the rational functions of ``\pi`` and at most seven parameters; ``\pi`` occupies one variable of the parameter field and lifts to the Symbolics constant. Its lifted output is not passed through `SQA.simplify`, which folds ``\pi`` into a float.
- Positive completion is explicit, never implicit in `floquet_expansion`, is defined only for Liouvillian expansions, and does not rewrite retained Floquet coefficients or micromotion.
- Kossakowski coordinates are relative to an ordered `JumpBasis`. Ordering is representation-significant even when two bases span the same subspace. The basis type and the exact extraction (`kossakowski`, `Lindbladian(S, basis)`, `support_basis`, `independent_basis`) belong to SecondQuantizedAlgebra (ADR 0014); FloquetExpansions adds methods only on `FloquetExpansion`.
- Raw expansions require an explicit `JumpBasis` for GKSL/Kossakowski coordinate extraction. Completed expansions store the finalized basis, so no-basis completed accessors are unambiguous.
- Automatic basis discovery is a symbolic convenience frontend whose output arity depends on runtime algebraic independence. The explicit-basis `positive_completion(expansion, algorithm, basis)` methods are the inference-oriented computational core.
- Automatic basis discovery starts from microscopic dissipative directions when provenance is available and appends the `support_basis` of every retained component, with independence taken modulo the identity by `independent_basis`. The microscopic directions are the Fourier harmonics of each channel operator, in SQA `channel_order_key` order of the channels and with the time average of each channel first, so a time-dependent channel never enters the basis with its phases.
- An empty `JumpBasis` is the basis of a coherent map. Automatic basis discovery returns one for an expansion without dissipative terms, so a coherent Liouvillian completes with no channels and a `0 × 0` Kossakowski matrix. GKSL extraction in the empty basis still rejects any two-sided residual.
- A completed result stores its finalized `JumpBasis` and caches the physical retained Kossakowski coefficients and coherent Hamiltonian used to construct its completed generator. A `JumpBasis` has no mutating API, so the stored basis cannot be changed from outside. Public accessors return defensive copies of matrices and channel vectors.
- Retained GKSL data in a fixed basis hold the coherent part already reattached to physical inverse-drive powers, while the dissipative coefficients stay order-resolved until completion finalization converts them to their physical series.
- `Gram()` is algebraic and must not require Hilbert-space/Liouville-space matrices, characteristic polynomials, symbolic eigendecomposition, or symbolic matrix square roots.
- `Spectral()` is a restricted perturbative spectral/HCM realization. The leading Kossakowski form must be diagonal in its basis, and retained corrections must already be diagonal within degenerate leading sectors. Automatic basis discovery does not diagonalize the leading Kossakowski form.
- `Lindbladian(cp)` and `channels(cp)` are defined only for completed and `GKSLNormalForm` expansions. The `SuperOperator` lowering of `Lindbladian(cp)`, with completeness expanded, equals `effective_generator(cp)`, and `channels(cp)` lists its `channel => rate` pairs in `channel_order_key` order.
- Algorithm-specific intermediate data are exposed through `factorization(cp)`; the normal completed physical interface is common to all completion algorithms.
- Expert API names marked `@public` are stable qualified interfaces, intentionally not widened into the ordinary export list.
- Dissipative quasienergy blocks use the energy-like Floquet-Liouville convention documented in ADR 0008; numerical vectorization remains an adapter concern.

Implementation-performance constraints such as solve-plan reuse and structured Hermitian elimination are owned by [`performance.md`](performance.md); the completion physics and representation contract are owned by ADR 0009.

When a proposed change conflicts with this map or an ADR, state the conflict and update or supersede the decision record before treating the new behavior as settled.

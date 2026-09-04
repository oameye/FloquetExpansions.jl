# CP-preserving completion of truncated Floquet Liouvillians

Finite-order Floquet Liouvillian expansions may leave the GKSL cone even when the microscopic driven model is Lindbladian. The package therefore supports an explicit positive-completion stage that selects a finite completely-positive continuation without changing any retained Van Vleck coefficient.

## Result model

Completion preserves the unified `FloquetExpansion` result model. A raw expansion carries an `Uncompleted` completion state. `positive_completion(vv, Gram())` or `positive_completion(vv, Spectral())` returns another `FloquetExpansion` carrying a concrete completed state. Completing an already completed expansion is an error; there is no `uncomplete` operation.

The retained `effective_components`, Hamiltonian components, gauge, and micromotion are not rewritten. The no-index `effective_generator(cp)` returns the finite completed generator, while `effective_component(cp, n)` continues to return the original retained order-`n` Floquet coefficient.

## Dissipative representation

Kossakowski coordinates are defined in an ordered `DissipativeFrame`. Ordering is semantically significant because coordinate matrices, deterministic pivot choices, and channel gauges depend on it even when two frames span the same operator subspace. Raw expansions require an explicit frame for Kossakowski extraction. Completion stores the finalized frame, so `dissipative_frame(cp)`, `kossakowski(cp)`, and `kossakowski_component(cp, n)` are unambiguous.

Automatic frame discovery starts from microscopic dissipative directions in user channel order and appends algebraically generated independent directions deterministically, removing dependence modulo the identity. `DissipativeFrame` stores the ordered operator data and coordinate information, not an eager `q × q` matrix of cross-dissipator Liouvillians.

## Microscopic provenance

`Liouvillian` and `PeriodicGenerator` remain basis-free algebraic representations. Arbitrary Liouvillian arithmetic does not propagate physical channel provenance.

The high-level `floquet_expansion(H, ωd, t, gauge, order; channels=...)` constructor may retain internal microscopic dissipative provenance before lowering the model to Fourier Liouvillian data. This provenance seeds automatic frame construction and sign assumptions during completion. Manually lowered `Liouvillian` or `PeriodicGenerator` inputs remain valid but receive no guaranteed microscopic provenance; unresolved signs are reported as completion conditions rather than guessed.

`jump(J, γ)` denotes a physical rate-weighted channel and therefore requires a real rate that is asserted nonnegative. A symbolic expression is asserted nonnegative as a whole. `collapse(L)` denotes a complete collapse operator and does not require extraction of a separate scalar rate.

## Completion algorithms

`Gram()` and `Spectral()` are public algorithm selectors with no required configuration fields in the initial API. A custom frame is supplied positionally:

```julia
positive_completion(vv, Gram())
positive_completion(vv, Spectral())
positive_completion(vv, Gram(), frame)
positive_completion(vv, Spectral(), frame)
```

The normal physical interface is common to both algorithms:

```julia
effective_generator(cp)
effective_component(cp, n)
hamiltonian(cp)
hamiltonian_component(cp, n)
kossakowski(cp)
kossakowski_component(cp, n)
dissipative_frame(cp)
channels(cp)
positivity_conditions(cp)
regularity_conditions(cp)
factorization(cp)
```

`channels(cp)` returns a channel representation that reconstructs the completed generator through `liouvillian(hamiltonian(cp); channels=channels(cp))`. Gram completion naturally returns complete collapse channels; spectral completion naturally returns rate-weighted jump channels. `channels` is intentionally undefined for an uncompleted expansion because the raw truncated effective generator need not admit a GKSL decomposition.

Algorithm-specific intermediate data are retained behind the single expert `factorization(cp)` accessor. Gram factorization data include the graded factor, channel onsets, and a compact active/dark diagnostic history. Spectral factorization data include completed branch rates, perturbative vectors/amplitudes, and onset information. The public physical API does not otherwise depend on the selected completion algorithm.

## Gram completion

For retained Kossakowski data

```math
d^{[N]}(ε)=\sum_{n=0}^{N} ε^n d^{(n)},
```

`Gram()` seeks a graded factor `B` satisfying

```math
\Pi_N(BB^\dagger)=d^{[N]},
```

and uses the untruncated finite Gram product `BB†` as the positive continuation. This guarantees positivity while changing only terms beyond the retained order.

The construction is algebraic. The leading Hermitian form is split into active and dark sectors by Hermitian congruence/elimination rather than eigendecomposition. The dark sector is the radical of the leading form and the active sector is the corresponding quotient. In reduced coordinates,

```math
d = \begin{pmatrix}A & X\\X^\dagger & C\end{pmatrix},
\qquad
\Sigma=C-X^\dagger A^{-1}X,
```

where the active block is factored by graded `LDL†` and scalar series square roots. The residual `Σ` is recursively classified by its first retained nonzero order. Negative retained directions are completion obstructions; positive even onsets open integer-power jump amplitudes; positive odd onsets signal Puiseux amplitudes. Residual orders beyond the truncation are closed because they are not constrained by the retained Floquet data.

`Gram()` must not require Hilbert-space or Liouville-space matrix representations, characteristic polynomials, symbolic eigendecomposition, or symbolic matrix square roots.

## Spectral completion

`Spectral()` provides the perturbative spectral/HCM realization in a restricted supplied frame. The leading Kossakowski form must be diagonal in that frame, and retained corrections must already be diagonal within degenerate leading sectors. Rayleigh–Schrödinger branch recursion is used only between distinct leading sectors. Completed nonnegative branch weights are reconstructed from square-completed rate series.

Automatic symbolic eigenframe discovery is not part of the initial implementation.

## Conditions and parameter strata

Completion distinguishes positivity conditions from regularity conditions. `positivity_conditions(cp)` reports unresolved inequalities required for nonnegative leading pivots or residual forms. `regularity_conditions(cp)` reports nonzero conditions needed for a fixed-rank symbolic factorization, for example denominators used in active-sector solves.

The first implementation uses deterministic structural sign reasoning and inherited physical rate assertions. A stronger SMT-backed sign oracle may be added later, but SymbolicSMT is not required by this design.

## Accuracy and scope

If completion is performed after retaining `N+1` Floquet coefficients, the finite completed generator agrees with the raw expansion through the same retained order. The retained micromotion is therefore reused unchanged at that accuracy.

Positive completion is nonunique beyond the retained order. Different frame orderings, congruence/pivot choices, Gram gauges, or spectral normalizations may produce different higher-order terms while preserving the same controlled perturbative data and complete positivity.

This construction does not imply convergence of the high-frequency expansion, uniqueness of the completion, or existence of an exact Floquet GKSL logarithm.

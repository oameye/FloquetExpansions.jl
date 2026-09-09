```@meta
CurrentModule = FloquetExpansions
CollapsedDocStrings = true
```

# Positive completion

A finite-order high-frequency expansion of a periodic Lindblad generator need not itself produce
a GKLS effective generator, even when the microscopic dynamics is Markovian and completely
positive [Schnell2021](@cite). Positive completion restores complete positivity without changing
the perturbative information fixed by the high-frequency expansion.

For a truncation through order ``N`` (`order = N + 1` in the API), write the retained Kossakowski
series as

```math
d^{[N]}(\epsilon)=\sum_{n=0}^{N}\epsilon^n d^{(n)},
\qquad \epsilon=\omega_d^{-1}.
```

A positive completion constructs a finite matrix ``\widetilde d^{[N]}(\epsilon)`` such that

```math
\widetilde d^{[N]}(\epsilon)\succeq0,
\qquad
\widetilde d^{[N]}(\epsilon)-d^{[N]}(\epsilon)=\mathcal O(\epsilon^{N+1}).
```

Thus the completion modifies only terms beyond the retained order. The perturbative effective
components remain unchanged, while [`effective_generator`](@ref) returns the completed finite
generator.

## Gram completion

[`Gram`](@ref) constructs a graded collapse-amplitude matrix

```math
B^{[N]}(\epsilon)=\sum_{p=0}^{N}\epsilon^p B^{(p)}
```

such that

```math
\Pi_N\!\left(B^{[N]}B^{[N]\dagger}\right)=d^{[N]}.
```

The finite completed Kossakowski matrix is the untruncated product

```math
\widetilde d^{[N]}=B^{[N]}B^{[N]\dagger}\succeq0.
```

The leading Hermitian form is split algebraically into active and dark sectors. The active block is
factorized with a graded ``LDL^\dagger`` construction and scalar series square roots; coupling to an
exactly dark sector is restored through linear solves and a Schur residual. No Kossakowski
eigendecomposition or Hilbert-space diagonalization enters this path.

When microscopic [`collapse`](@ref) or [`jump`](@ref) channels are available, automatic frame
construction preserves their order and only appends additional generated directions needed by the
retained Floquet components. An explicit [`DissipativeFrame`](@ref) fixes the representation.

```julia
cp = positive_completion(vv, Gram())
# or
cp = positive_completion(vv, Gram(), frame)

kossakowski(cp)
dissipative_frame(cp)
channels(cp)
factorization(cp)
```

[`factorization`](@ref) returns a [`GramFactorization`](@ref). Its `amplitudes` contain the physical,
drive-frequency-reattached amplitude grades. `positivity_conditions` records scalar assumptions of
the form ``p\ge0``; `regularity_conditions` separately records nonzero conditions needed to remain
on the selected fixed-rank symbolic stratum.

The first Gram stage handles a regular leading-rank sector and exactly dark residual directions.
A new dissipative channel that opens only at a higher perturbative grade is classified before the
recursive continuation step; that recursive onset filtration is a separate completion stage.

## Spectral completion

[`Spectral`](@ref) works in a restricted dissipative frame in which the leading Kossakowski form is
diagonal. It follows perturbative decay-rate branches and completes each retained branch rate by an
HCM square construction. Mixing inside an unresolved degenerate leading sector is rejected rather
than resolved by an implicit symbolic diagonalization.

Odd branch-rate onsets remain valid rate-weighted jumps. [`SpectralFactorization`](@ref) records the
branch rates, normalized perturbative vectors, onset orders, and whether a rate-folded collapse
amplitude would require a Puiseux onset.

## Micromotion

The retained kick operator does not need to be modified. If
``\mathcal L_{\rm CP}^{[N]}-\mathcal L_{\rm eff}^{[N]}=\mathcal O(\epsilon^{N+1})``, then the same
retained micromotion gives a completed approximation with the same accuracy,

```math
\Phi_{\rm CP}^{[N]}(t)
=
e^{\mathcal K^{[N]}(t)}
e^{t\mathcal L_{\rm CP}^{[N]}}
e^{-\mathcal K^{[N]}(0)}.
```

A change of micromotion is only required if one asks for exact equality with the original finite
truncated propagator rather than agreement through the retained order.

## API

```@docs
Completion
Uncompleted
CompletionAlgorithm
Gram
Spectral
CompletionFactorization
CompletionObstruction
FractionalJumpOnset
GramFactorization
SpectralFactorization
positive_completion
dissipative_frame
channels
positivity_conditions
regularity_conditions
factorization
```
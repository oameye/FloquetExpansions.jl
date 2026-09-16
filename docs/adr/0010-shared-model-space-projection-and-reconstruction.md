# Shared model-space projection and reconstruction boundaries

`FloquetExpansions.jl` treats model-space projection as a reusable internal construction that is distinct from the physical/output reconstruction built on top of it. This decision records the package-level conventions settled in #109 and #110.

The package architecture must keep the following axes separate:

```text
model-space choice
    ≠ projection core / complement solve
    ≠ expansion algorithm
    ≠ normalization / reconstruction semantics
    ≠ gauge / physical representative
    ≠ transformation representation
    ≠ algebra specialization
```

## Shared projection vocabulary

For a retained model space `P` and complement `Q = 1 - P`, use

```math
\Omega = P + Y,
\qquad
Y = Q\Omega P,
```

with the model-space embedding / Bloch wave operator defined by

```math
\mathbb G\Omega = \Omega B,
\qquad
P\Omega P = P.
```

`Ω` is the model-space embedding or Bloch wave operator. `Y` is the invariant graph / complement dressing. `B` is the intermediate-normalized model-space generator. A **complement solve** is the generic operation used to determine the `Q`-space response; call it a reduced resolvent only when it is genuinely resolvent-like.

Bloch wave operator, invariant graph, Feshbach/Bloch-Horowitz elimination, Brillouin-Wigner, and Van Vleck are related constructions but are not interchangeable names. Van Vleck is an output normalization/representative, not the shared projection object itself.

For ordinary Floquet HFE,

```text
P = retained/slow Floquet block
Q = virtual Floquet complement
Q-solve = inverse fast harmonic operator / Floquet reduced resolvent
```

For higher-order QHB,

```text
P = explicitly retained carrier model space
Q = generated omitted harmonics
Q-solve = physical omitted-harmonic response solve
```

The same projection vocabulary may be reused, but the physical model-space choice, complement solve, small parameter, and complexity need not be the same. In particular, the `O(N²)` perturbative-order count established for the nonresonant Floquet recurrence must not be transferred automatically to QHB.

## Internal data flow

The intended internal flow is

```text
canonical physical/process data
        ↓
model-space problem (P/Q + complement solve + perturbative bookkeeping)
        ↓
projection solution (Ω/Y, B, provenance/diagnostics)
        ↓
output-specific reconstruction
```

Physical-process provenance, harmonic/grade kinematics, `P/Q` state, and output-reconstruction metadata are distinct axes. Where physical provenance is available, lowering to projection data must not erase the identity of a microscopic channel, physical parent process, adjoint partner, or generated descendant.

This is an internal architectural contract, not a requirement to expose public `PhysicalProcess`, `ProjectionProblem`, or `ProjectionSolution` types. Public abstractions are introduced only when a validated implementation requires them.

## Canonical Van Vleck reconstruction

The canonical Van Vleck branch consumes shared projection data through

```text
shared projection data
    ↓
model-space canonical normalization
    ↓
periodic transformation normalization
    ↓
connected / noncommutative logarithm
    ↓
static frame/similarity correction
    ↓
package VanVleck() effective generator + zero-average single kick
```

Use **canonical normalization** for the Okubo/des-Cloizeaux/model-space normalization step. Use **static frame/similarity correction** for the additional static transformation required to reach the package representative. Reserve **Van Vleck gauge** for the final zero-average micromotion condition.

For generic Liouvillians the model-space transformation is a similarity, not a unitary transformation and not generally a completely-positive map. Exact package `VanVleck()` equivalence therefore requires both the effective generator and the retained zero-average micromotion coefficients, not spectral agreement alone.

## CP-preserving HFE reconstruction

The built-in CP-preserving HFE is a sibling reconstruction with different finite-order semantics:

```text
physical Hamiltonian + jump/amplitude data
    ↓
shared coherent/projection solve where appropriate
    ↓
physical jump/amplitude transport
    ↓
finite Gram/Kraus reconstruction
    ↓
manifestly GKSL/CP effective representation
```

The construction approximates physical amplitudes before forming the positive quadratic object. For a retained amplitude factor `A^[N]`, the finite output uses

```math
\widetilde d^{[N]} = A^{[N]} A^{[N]\dagger} \succeq 0
```

without truncating the final square again. The resulting higher-order tail is part of the finite CP approximant.

The CP branch must not compute a generic truncated Liouvillian Van Vleck expansion and then repair or complete it. It must not be represented as a `Gauge`, and it must not call `positive_completion` internally. The terms **CP amplitude reconstruction** and **Gram/Kraus reconstruction** are used for this native HFE step.

The one-dissipator static relation to canonical Liouvillian Van Vleck is described by a **static model-space/frame similarity correction**. That similarity is not generally CP and therefore does not turn canonical Liouvillian Van Vleck into the physical CP representation.

## Relationship to positive completion

The positive-completion architecture in ADR 0009 remains a separate package capability.

```text
positive completion
    already-truncated Van Vleck Kossakowski data
        → positive higher-order continuation preserving retained coefficients

native CP HFE
    physical Hamiltonian + amplitudes/jumps
        → amplitude-level HFE
        → finite Gram/Kraus reconstruction
```

The two paths may reuse representation-agnostic infrastructure such as exact series algebra, dissipative coordinates, Gram-factor storage, and channel reconstruction, but they do not share defining semantics.

## API and implementation constraints

Existing `floquet_expansion(..., VanVleck(), order)` behavior remains unchanged unless a separately validated backend is deliberately selected later.

Do not add public algorithm, reconstruction, or transformation hierarchies solely for future-proofing. When #97, #105, or another concrete implementation requires a seam, add the smallest static-dispatch abstraction that the implementation actually exercises.

Hot-path data must remain concrete and inference-friendly. Avoid `Any`, abstract-typed storage, runtime string/symbol dispatch, and dense materialization when a compact natural representation is sufficient. Generic formal-series, harmonic-domain, projection-plan, workspace, and cache machinery should be reused across output branches where the mathematics is genuinely shared; normalization- or reconstruction-specific semantics must remain in their own layers.

Performance accounting must separate shared projection cost from output reconstruction cost. Recurrence-level operation counts do not imply total symbolic/runtime complexity.

## Consequences

#97 may build canonical normalization/logarithmic reconstruction on the shared projection solution without defining the CP path. #105 may reuse the coherent/projection and representation-agnostic series kernels without passing through generic Liouvillian Van Vleck. Higher-order QHB may reuse the `P/Q` projection machinery with a different model space and complement solve.

#98 should extract a genuinely common compiled projection IR from these concrete consumers rather than force all outputs into one representation in advance. #72 owns any later public/static implementation seam once a real second backend or reconstruction requires it.

## Gate

Until the concrete consumers land, this ADR is held by the architecture/specification decisions in #109 and #110. Implementations under #97, #105, #113/#115, #98, and #72 must preserve these distinctions and add exact behavior, compiler, and performance regressions at the layer they introduce.

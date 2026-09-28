# FloquetExpansions.jl

Shared vocabulary for symbolic high-frequency expansions of periodically driven quantum systems.

The authoritative research/dependency roadmap is #45. The package currently has a production ordinary-HFE + post-hoc-completion stack and two distinct open-system research branches:

```text
ordinary/generic HFE
    |
    +-- post-hoc positive completion          #54/#313
    |
    +-- native GKLS HFE                       #350, Phase I
    |
    `-- Floquet–Kraus / CPTP period map       #168/#169, Phase II
```

The phrase **native CP HFE** is deprecated because it conflates the last two branches.

## Dynamics

**Periodically driven Lindblad system**:
An open quantum system whose Lindblad generator varies periodically with time.

**Lindblad / GKLS generator**:
A density-operator generator written as a Hermitian Hamiltonian action plus dissipative channels with a positive Kossakowski form.
_Avoid_: using this term for a finite-order effective result unless its GKLS form is known.

**Liouvillian**:
A linear map that generates or approximates density-operator evolution; it need not have GKLS form.
_Avoid_: treating every Liouvillian as completely positive or Markovian-embeddable.

**Hamiltonian action**:
The coherent contribution ``-i[H, ρ]`` to a density-operator generator.

**Dissipator**:
The map ``D[L](ρ) = LρL† - (L†Lρ + ρL†L)/2`` associated with a collapse operator.

**Collapse operator**:
An operator used directly in ``D[L]``; any amplitude or phase belonging to that channel is included in ``L``.

**Rate-weighted jump channel**:
A bare jump operator and a separate real, nonnegative scalar rate contributing ``γD[J]``; the rate is not part of ``J``. A symbolic rate expression is asserted nonnegative as a whole.

**Time-dependent dissipative channel**:
A dissipative channel whose operator or scalar rate varies periodically with the drive.

## Floquet expansion

**Periodic generator**:
A time-periodic generator expressed through Fourier harmonics in a drive-frequency basis.
_Avoid_: `PeriodicOperator`, the superseded Hamiltonian-only name.

**Van Vleck expansion**:
The high-frequency expansion in the van Vleck gauge, separating a periodic generator into a time-independent effective generator and periodic micromotion.
_Avoid_: treating van Vleck as a computational algorithm; it is a Floquet representative/gauge.

**Van Vleck pair**:
The effective generator and micromotion generator of the van Vleck gauge, taken together. The zero-average condition fixes the pair, which is why independent expansion algorithms for that gauge must return the same retained representative.

**Floquet gauge / representative**:
The convention fixing the static freedom in the periodic/effective decomposition. The van Vleck representative uses zero-average logarithmic micromotion; a Floquet–Magnus/stroboscopic representative instead fixes an endpoint/initial phase.
_Avoid_: conflating Floquet gauge with the native-GKLS static normal-form freedom inside a fixed Floquet branch.

**Expansion algorithm**:
The procedure that computes the coefficients of a chosen Floquet representative. The package algorithm axis is

```text
HoriDeprit()
BlochFeshbach()
```

Algorithms for the same representative differ in intermediate objects and cost, not in the retained physical representative after normalization.
_Avoid_: calling an algorithm a gauge, a positivity policy, or a physical reconstruction.

**Hori–Deprit algorithm**:
The Lie-transform expansion algorithm, solving the homological equations order by order through nested commutators.

**Bloch/Feshbach algorithm**:
The wave-operator/invariant-graph algorithm. It solves the model-space problem in intermediate normalization and then converts to the requested canonical Floquet representative.

**Model space**:
The retained subspace selected by a projector ``P``. Its complement is ``Q=1-P``. For ordinary periodic HFE, ``P`` is the retained/zero Floquet block. Other consumers such as QHB choose a different physical model space.

**Floquet homological inverse**:
The complement solve that eliminates nonzero Floquet harmonics. For the periodic Fourier convention ``e^{-imτ}``, schematically

```math
∂_τ^{-1}e^{-imτ} = \frac{i}{m}e^{-imτ},\qquad m\ne0,
```

or equivalently the corresponding Sambe/Feshbach ``Q`` resolvent.

**Wave operator / invariant graph**:
A near-identity periodic map carrying model-space solutions into the full driven problem. In Bloch/Feshbach intermediate normalization it is the natural object before canonical Floquet normalization.

**Bloch effective generator**:
The constant model-space generator paired with an intermediate-normalized Bloch wave operator. It is similar to the canonical van Vleck effective generator but is a different representative.
_Avoid_: calling it the final package effective generator before normalization.

**Static factor**:
The constant factor ``N`` that renormalizes a Bloch wave operator to the van Vleck representative through ``ΩN=e^Λ`` with ``⟨Λ⟩=0``. The canonical effective generator is the corresponding similarity transform ``N⁻¹BN``.

**Connected logarithm**:
The periodic logarithm ``Λ`` of the normalized transformation ``ΩN``. In the canonical branch it is organized as a connected/Lie series.
_Avoid_: calling the raw logarithm of ``Ω`` itself the canonical micromotion logarithm.

**Floquet expansion**:
The finite-order result of applying a high-frequency expansion to a periodic generator, including retained effective-generator coefficients, micromotion, and completion state where applicable.

**Effective generator**:
The time-independent generator used in the chosen effective description. Its semantics depend on the branch:

- ordinary/raw HFE: algebraic truncation of the chosen Floquet representative;
- post-hoc completion: selected finite positive continuation preserving all retained raw coefficients;
- native GKLS HFE: the static GKLS representative solved together with micromotion inside the recurrence.

**Effective component**:
A retained order-by-order coefficient of the effective generator. Post-hoc positive completion does not rewrite retained raw components. Native GKLS HFE is different: its retained coefficient can move within the static Floquet similarity orbit because the compensating change is carried by micromotion.

**Finite-order effective generator**:
A truncated algebraic result. For a generic Liouvillian HFE it need not retain GKLS form or complete positivity.

**Micromotion**:
The periodic transformation relating the original time-dependent dynamics to an effective static description. For generic dissipative systems it is a similarity map, not necessarily a CP map and not a unitary kick.

**Dissipative quasienergy**:
A generally complex eigenvalue of the energy-like Floquet-Liouville operator; its real part describes oscillation while its imaginary part describes decay or growth, modulo the drive frequency.

## Open-system constructions

These three constructions are deliberately distinct.

### Post-hoc positive completion — #54/#313

**Positive completion**:
An explicit post-processing step applied to an already-computed raw canonical HFE. It selects a finite positive Kossakowski continuation while preserving every retained raw Floquet coefficient and the retained micromotion. It is generally nonunique beyond retained order.

```text
raw retained canonical HFE
    -> positive higher-order continuation
    -> retained raw coefficients unchanged
```

It is not the native GKLS HFE and it is not a Floquet–Kraus/CPTP map construction.

### Native GKLS HFE — #350, Phase I

**Native GKLS HFE**:
The generator-native high-frequency branch that solves the retained ``m=0`` Floquet coefficient directly in GKLS coordinates while the ordinary Floquet homological recurrence is being advanced.

At order ``n`` the static dark equation is schematically

```math
Δ_n + Φ_n(S_n) = P_n,\qquad P_n\succeq0,
```

with

```math
Φ_n(S)=(π_n\otimes\bar π_n)c([\bar{\mathcal L},S]).
```

The unknowns include the intrinsic static normal-form coefficient, effective Hamiltonian, corrections to active jump amplitudes, and newly born channels when required.

**Static GKLS homological operator**:
``Φ_n`` (or its generalized inverse ``Φ_n^#``) acting inside the already-retained Floquet block. It fixes residual static similarity so that the effective coefficient admits GKLS coordinates.
_Avoid_: calling this a second Floquet ``P/Q`` projection. The primary Floquet homological inverse and this retained-block static solve are different operations.

**Intrinsic static GKLS normal-form freedom**:
The residual static similarity used by #350 after the Floquet model space has been selected. It is defined from the total kick/static similarity, not from an algorithm-specific HD/BF slot.
_Avoid_: calling it simply “the Floquet gauge”.

**Graded jump amplitude / Gram column**:
A generator coordinate carried by #350, possibly with half-integer amplitude onset while the assembled generator has an ordinary integer-power expansion.
_Avoid_: calling these Kraus operators of the one-period dynamical map.

The Hamiltonian-only limit is a hard reduction requirement: #350 must reduce coefficient by coefficient to ordinary Hamiltonian Hori–Deprit/Bloch–Feshbach/van-Vleck HFE.

### Floquet–Kraus / CPTP propagator — #168/#169, Phase II

**Floquet–Kraus propagator**:
The map-first open-system construction in which physical one-branch/Stinespring amplitudes are Floquet-expanded/projected before left/right pasting and trace-preserving reconstruction. Its primary output is a CPTP one-period map.

```text
physical one-branch amplitudes
    -> Floquet/Sambe amplitude-space P/Q problem
    -> irreducible half-chain kernels
    -> TP-coupled physical reconstruction
    -> CPTP period map
```

**Kraus / Stinespring amplitude**:
A one-branch dynamical-map amplitude used by the Phase-II construction.
_Avoid_: identifying it with a #350 static jump/Gram column.

**CPTP map**:
A completely positive, trace-preserving dynamical map. CPTP is a map-level structural guarantee, not a statement that its logarithm is GKLS.

**GKLS embeddability**:
The additional property that a CPTP period map, in an appropriate Floquet branch/representative, can be written as the exponential of a time-independent GKLS generator. This is not universal.

When comparing Phase II to #350, the raw connected logarithm of a stroboscopic period map must first be put in the same Floquet representative. In particular, endpoint micromotion/static-gauge conversion is required before comparing a stroboscopic/Floquet–Magnus logarithm with a canonical van-Vleck/native-GKLS representative.

## CP, TP, and GKLS

Use these terms separately:

**CP**: complete positivity of a map/operator-sum construction.

**TP**: trace-preserving normalization/completeness.

**CPTP**: both CP and TP; a physical quantum channel/map guarantee.

**GKLS**: a time-local Markovian generator representation with Hermitian Hamiltonian and positive dissipative form.

_Avoid_: using CP, CPTP, and GKLS interchangeably.

## GKSL coordinates and positive completion

**Dissipative frame**:
An ordered finite set of operator directions ``(F₁, …, F_q)`` used to represent the dissipative Hermitian form. The frame need not be complete, traceless, or Hilbert-Schmidt orthonormal. Ordering is part of the representation and affects coordinate matrices, deterministic pivot choices, and channel gauge.

**Kossakowski matrix**:
The Hermitian coordinate matrix ``d`` of the dissipative form in a specified dissipative frame. It is representation dependent, while positive semidefiniteness and inertia are preserved by nonsingular congruence.

**Completion state**:
The type-level state carried by a `FloquetExpansion`, distinguishing an uncompleted expansion from a selected post-hoc positive continuation. Completing an already completed expansion is an error.

**Gram completion**:
An algebraic positive-completion algorithm constructing a graded Gram factor ``B`` with retained matching ``Π_N(BB†)=d^[N]`` and using the untruncated finite product ``BB†`` as the positive continuation. It uses Hermitian congruence, graded ``LDL†``, and recursive Feshbach/Schur reduction rather than symbolic eigendecomposition.

**Spectral completion**:
A perturbative spectral/HCM positive completion following decay-rate branches in a restricted spectral frame and square-completing retained rates. It is an independent realization/validation oracle rather than a prerequisite for Gram completion.

**Gram factor**:
A matrix ``B`` satisfying ``d=BB†`` for a positive Kossakowski form. Right-unitary/isometric rotations of its columns are jump-channel gauge transformations and leave ``d`` unchanged.

**Active dissipative sector**:
The nondegenerate quotient of the leading dissipative Hermitian form after removing its radical. It contains channels already open at the current perturbative grade.

**Dark sector**:
The radical of the current leading dissipative Hermitian form. Its directions have no rate at that grade and are resolved at higher grades through the reduced residual form.

**Feshbach residual**:
For a reduced block form ``d=[A X; X† C]`` with active ``A``, the induced dark-sector Hermitian form ``Σ=C-X†A⁻¹X``. It determines whether dark directions open, remain unresolved beyond truncation, or obstruct a positive continuation.

**Dissipative onset filtration**:
The nested sequence of dark sectors exposed by recursively resolving Feshbach residuals. Its successive quotients collect channels that open at successive perturbative orders.

**Positivity condition**:
A symbolic inequality required for a completion to be nonnegative when it cannot be established structurally or from inherited physical-rate assertions.

**Regularity condition**:
A nonzero condition defining the fixed-rank parameter stratum on which a symbolic factorization is valid.

**Microscopic dissipative provenance**:
Internal information retained only by the high-level physical `floquet_expansion(...; channels=...)` construction, used to seed automatic dissipative-frame discovery and physical rate assumptions. Generic `Liouvillian` and `PeriodicGenerator` algebra remains provenance-free.

**Completion factorization**:
Algorithm-specific diagnostic data retained by a post-hoc completed Floquet expansion and exposed through `factorization`.

## Gauge freedoms that must not be conflated

1. **Floquet representative/static gauge** — van Vleck, Floquet–Magnus/stroboscopic, etc.
2. **Native-GKLS static normal-form freedom** — #350 retained-block static similarity.
3. **Jump-channel unitary gauge** — rotations among static GKLS jump/Gram columns.
4. **Kraus/channel unitary gauge** — equivalent Kraus representations of a CPTP map.
5. **CK/purification amplitude gauge** — amplitude-space representation freedom in the Phase-II construction.
6. **Affine Lindblad representation freedom** — shifts of jump operators by identity with Hamiltonian compensation where applicable.

## Roadmap phase names

```text
Phase I  = #350 native GKLS HFE
Phase II = #168/#169 Floquet–Kraus / CPTP propagator
```

The historical #104/#105/#122 coherent-frame work is the **one-dissipator/coherent-frame oracle**, not “Phase I”. Post-hoc completion is an independent package capability, not a roadmap phase.

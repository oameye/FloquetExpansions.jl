# FloquetExpansions.jl

Shared vocabulary for symbolic high-frequency expansions of periodically driven quantum systems.

## Dynamics

**Periodically driven Lindblad system**:
An open quantum system whose Lindblad generator varies periodically with time.

**Lindblad generator**:
A density-operator generator written as a Hamiltonian action plus dissipative channels with physically valid rates.
_Avoid_: using this term for a finite-order effective result unless its Lindblad form is known.

**Liouvillian**:
A linear map that generates or approximates density-operator evolution; it need not have Lindblad form.
_Avoid_: treating every Liouvillian as completely positive.

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
_Avoid_: treating van Vleck as a method applied in some other gauge; it is the gauge, and the method is the expansion algorithm.

**Van Vleck pair**:
The effective generator and the micromotion generator of the van Vleck gauge, taken together. The zero-average condition fixes the pair uniquely, which is why every expansion algorithm for that gauge must return it.

**Gauge**:
The convention fixing the free integration constant of the micromotion generator, and with it which effective generator and micromotion a Floquet expansion represents. The van Vleck gauge requires the micromotion generator to have zero period average.

**Expansion algorithm**:
The procedure that computes the coefficients a gauge fixes. Algorithms for the same gauge produce identical retained effective components and micromotion; they differ only in intermediate quantities and cost.
_Avoid_: calling an algorithm a gauge, or confusing it with a positive-completion algorithm, which changes the finite effective generator.

**Hori–Deprit algorithm**:
The Lie-transform expansion algorithm, which solves for the micromotion generator order by order through nested commutators of drive harmonics.

**Bloch/Feshbach algorithm**:
The wave-operator expansion algorithm, which solves the Bloch equation for a periodic wave operator and a Bloch effective generator, then normalizes both to the gauge's micromotion and effective generator.

**Wave operator**:
A periodic map with unit period average that carries solutions of a constant generator into solutions of the driven dynamics.

**Bloch effective generator**:
The constant generator paired with a wave operator. It is similar to the van Vleck effective generator and shares its quasienergies, but it is a different representative and is not Hermitian in general for a Hamiltonian.
_Avoid_: calling it the effective generator, which is always the gauge's representative.

**Static factor**:
The constant ``N`` that renormalizes a wave operator to the van Vleck gauge through ``ΩN = e^{Λ}`` with ``⟨Λ⟩ = 0``. The van Vleck effective generator is its similarity transform ``N⁻¹BN`` of the Bloch effective generator.

**Connected logarithm**:
The periodic logarithm ``Λ`` of the renormalized wave operator ``ΩN``. It has zero period average and contains only nested commutators of drive harmonics, so it is the van Vleck micromotion generator.
_Avoid_: calling it the logarithm of the wave operator, which differs from it by the static factor.

**Floquet expansion**:
The finite-order result of applying a high-frequency expansion to a periodic generator, including its retained effective-generator coefficients, micromotion, and completion state.

**Effective generator**:
The time-independent generator that approximates the slow dynamics after the expansion. For an uncompleted Floquet expansion this is the raw algebraic truncation; for a positively completed expansion it is the selected finite positive continuation.

**Effective component**:
A retained order-by-order coefficient of the Floquet effective generator. Positive completion does not rewrite retained effective components.

**Finite-order effective generator**:
A truncated algebraic result that may not retain Lindblad form or complete positivity.

**Micromotion**:
The periodic transformation that relates the original time-dependent dynamics to the effective dynamics; its generator need not be unitary for dissipative systems. Positive completion of the effective generator does not modify retained micromotion through the controlled perturbative order.
_Avoid_: calling dissipative micromotion a unitary kick operator.

**Dissipative quasienergy**:
A generally complex eigenvalue of the energy-like Floquet-Liouville operator for a periodically driven open system; its real part describes oscillation while its imaginary part describes decay or growth, and it is defined modulo the drive frequency.
_Avoid_: assuming dissipative quasienergies are real or that every Liouvillian is completely positive.

## GKSL coordinates and positive completion

**Dissipative frame**:
An ordered finite set of operator directions ``(F₁, …, F_q)`` used to represent the dissipative Hermitian form. The frame need not be complete, traceless, or Hilbert-Schmidt orthonormal. Ordering is part of the representation and affects coordinate matrices, deterministic pivot choices, and channel gauge.
_Avoid_: treating two differently ordered frames as identical merely because they span the same subspace.

**Kossakowski matrix**:
The Hermitian coordinate matrix ``d`` of the dissipative form in a specified dissipative frame. It is representation dependent, while positive semidefiniteness and inertia are preserved by nonsingular congruence.

**Positive completion**:
An explicit post-processing step that selects a finite positive Kossakowski continuation while preserving every retained Floquet coefficient. It is not implicit in `floquet_expansion` and is generally nonunique beyond the retained order.

**Completion state**:
The type-level state carried by a `FloquetExpansion`, distinguishing an uncompleted expansion from a selected positive continuation. Completing an already completed expansion is an error.

**Gram completion**:
An algebraic positive-completion algorithm that constructs a graded Gram factor ``B`` with retained coefficient matching ``Π_N(BB†) = d^[N]`` and uses the untruncated finite product ``BB†`` as the positive continuation. The construction uses Hermitian congruence, graded ``LDL†`` factorization, and recursive Feshbach/Schur reduction rather than symbolic eigendecomposition.

**Spectral completion**:
A perturbative spectral/HCM positive completion that follows decay-rate branches in a restricted spectral frame and square-completes retained rates. It is an independent realization and validation oracle rather than a prerequisite for Gram completion.

**Gram factor**:
A matrix ``B`` satisfying ``d = BB†`` for a positive Kossakowski form. Right-unitary/isometric rotations of its columns are jump-channel gauge transformations and leave ``d`` unchanged.

**Active dissipative sector**:
The nondegenerate quotient of the leading dissipative Hermitian form after removing its radical. It contains channels already open at the current perturbative grade.

**Dark sector**:
The radical of the current leading dissipative Hermitian form. Its directions have no rate at that grade and are resolved at higher grades through the reduced residual form.

**Feshbach residual**:
For a reduced block form ``d = [A X; X† C]`` with active ``A``, the induced dark-sector Hermitian form ``Σ = C - X†A⁻¹X``. It determines whether dark directions open, remain unresolved beyond truncation, or obstruct a positive continuation.

**Dissipative onset filtration**:
The nested sequence of dark sectors exposed by recursively resolving Feshbach residuals. Its successive quotients collect channels that open at successive perturbative orders and replace perturbative eigenvalue branches as the general algebraic organizing structure.

**Positivity condition**:
A symbolic inequality required for a completion to be nonnegative when it cannot be established from structural algebra or inherited physical rate assertions.

**Regularity condition**:
A nonzero condition defining the fixed-rank parameter stratum on which a symbolic factorization is valid, for example an active pivot used in division or a scalar series square root.

**Microscopic dissipative provenance**:
Internal information retained only by the high-level physical `floquet_expansion(...; channels=...)` construction, used to seed automatic dissipative-frame discovery and physical rate assumptions. Generic `Liouvillian` and `PeriodicGenerator` algebra remains provenance-free.

**Completion factorization**:
Algorithm-specific diagnostic data retained by a completed Floquet expansion and exposed through `factorization`. Gram data include the graded factor, onset information, and compact active/dark history; spectral data include completed branch rates, perturbative vectors/amplitudes, and onset information.

## Harmonic balance

**Quantum harmonic balance**:
A description of one driven resonator through several carrier frequencies at once. The carrier modes are bookkeeping for one physical mode, not additional microscopic resonators, and they do not carry their own baths.
_Avoid_: QHB as a gauge or an expansion algorithm; it is a separate model construction.

**Carrier**:
A frequency ``ω_j`` at which the physical mode is resolved. Carrier frequencies are symbolic expressions with exact coefficients, and relations between carriers are written into those expressions.

**Carrier mode**:
One of the two ladder modes of a carrier: the co-rotating mode ``a_{j+}``, which enters with the phase ``e^{-iω_j t}``, and the counter-rotating mode ``a_{j-}``, which enters with ``e^{+iω_j t}``.

**Carrier embedding**:
The map ``I_t(a) = λ_H^{-1/2} \sum_j (e^{-iω_j t} a_{j+} + e^{+iω_j t} a_{j-})`` from the physical mode into the carrier modes. It is an algebra homomorphism, a unitary change of frame on the physical mode extended by dark partners.

**Harmonic normalization**:
The number ``λ_H = 2N_h`` of carrier modes for ``N_h`` carriers. It fixes the embedding weight ``λ_H^{-1/2}`` and multiplies the resonant part of the carrier generator.

**Frame Hamiltonian**:
The rotating-frame gauge term ``H_{fr} = -\sum_j ω_j (n_{j+} - n_{j-})`` of the carrier generator.

**Resonant projection**:
The projection ``P`` that keeps the processes whose total frequency expands to zero. For a Liouvillian it acts on each ``ρ ↦ AρB`` term as a whole.

**Unreduced harmonic balance**:
The carrier generator ``λ_H P[I_t(\cdot)] + H_{fr}`` with both carrier modes of every carrier kept. It is exact for quadratic Hamiltonians and linear dissipators and approximate for nonlinear terms.
_Avoid_: QHB-RWA, the later reduction that keeps one mode per carrier.

**Reconstruction**:
Reading a physical observable from a carrier state through ``⟨A⟩(t) = \mathrm{Tr}[I_t(A)ρ]``.
_Avoid_: micromotion, which belongs to a Floquet expansion.

**Bright and dark sectors**:
The combination of carrier modes that the embedding assigns to the physical mode at time ``t`` is bright. Every carrier combination the reconstruction does not see is dark.

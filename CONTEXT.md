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
A high-frequency expansion that separates a periodic generator into a time-independent effective generator and periodic micromotion in a chosen gauge.

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

## Model-space projection and reconstruction

**Model space**:
The retained sector selected by a projector ``P`` in a projection/invariant-subspace construction. Its physical meaning depends on the method: for ordinary Floquet HFE it is the retained/slow Floquet block; for higher-order QHB it is the explicitly retained carrier space.
_Avoid_: assuming every model space has the same small parameter, complement solve, or complexity.

**Complement**:
The eliminated or response sector ``Q = 1 - P`` associated with a chosen model space.

**Model-space embedding / Bloch wave operator**:
The intermediate-normalized map ``Ω = P + Y`` satisfying ``𝔾Ω = ΩB`` and ``PΩP = P``. It embeds retained model-space data into the dressed full representation.
_Avoid_: using `Van Vleck`, `Brillouin-Wigner`, or `Feshbach` as interchangeable names for ``Ω``.

**Invariant graph**:
The complement dressing ``Y = QΩP`` of the model-space embedding.

**Intermediate-normalized model-space generator**:
The retained generator ``B`` appearing in ``𝔾Ω = ΩB`` before output-specific canonical normalization or physical reconstruction.
_Avoid_: calling ``B`` the package `VanVleck()` generator without the required normalization/gauge conversion.

**Complement solve**:
The method-specific operation that determines the ``Q``-space response. It may be a fast-harmonic inverse, a reduced resolvent, a Sylvester solve, or an omitted-harmonic response solve.
_Avoid_: calling every complement solve a resolvent when no resolvent is actually formed.

**Reconstruction**:
The generic operation that converts model-space/projection data into a physical or package output representation. `Micromotion`/`kick` is reserved for Floquet representatives that actually use that structure; QHB readout is physical reconstruction, not Floquet micromotion.

**Observable dressing**:
The projection-induced reconstruction of physical observables or readout from retained model-space data. It is the generic term for this operation across Floquet HFE and higher-order QHB.
_Avoid_: calling generic observable dressing `micromotion` when the construction is not a Floquet representative.

**Canonical normalization**:
The Okubo/des-Cloizeaux/model-space normalization that converts an intermediate Bloch/Feshbach representation toward the canonical Van Vleck representative.

**Static frame/similarity correction**:
A static model-space transformation used to relate normalized representatives or impose the final package representative. For a generic Liouvillian it is a similarity, not generally a unitary or completely-positive map.

**Van Vleck gauge**:
The final zero-average micromotion condition of the package Van Vleck representative.
_Avoid_: using `gauge` as a synonym for every model-space normalization or reconstruction step.

**CP-preserving high-frequency expansion**:
A native high-frequency construction that transports physical jump/amplitude data and forms a finite positive Gram/Kraus object, rather than completing an already-truncated Liouvillian Van Vleck generator.
_Avoid_: treating CP preservation as a `Gauge` or as implicit positive completion.

**CP amplitude reconstruction**:
The output reconstruction that forms the finite positive object from retained physical amplitudes or equivalent Kraus/Gram factors. The final positive square is not re-truncated merely to recover a polynomial Kossakowski series.

**Physical-process provenance**:
Internal identity and origin data associated with a physical process before output-specific lowering, such as its operator/map payload, coefficient, harmonic/grade, adjoint partner, microscopic channel/bath identity where present, and parent/descendant relation when generated. Provenance is distinct from harmonic kinematics, ``P/Q`` state, and reconstruction metadata.

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

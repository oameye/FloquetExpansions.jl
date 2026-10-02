# The GKSL normal form solves its jump amplitudes inside the BF/HD recurrence

The generator-native GKLS expansion does not compute an ordinary Liouvillian Van Vleck expansion and repair its Kossakowski matrix afterwards. At every order the static slot of the Bloch–Feshbach or Hori–Deprit recurrence is chosen so that the order-``n`` coefficient already has the graded Gram form ``c(E_n) = K_n^< + K_n^{\rm lin} + P_n``, and the accepted ``E_n`` is fed back into the same recurrence. The static gauge ``S_n``, the coherent part ``H_n``, and the graded jump amplitudes are therefore state of the recurrence. Every finite truncation is GKLS by construction, and positive completion (ADR 0009) is never called.

The implementation lives in `src/native/`, after `gksl/` and independent of `completion/`:

- `GradedChannel{T}` holds an amplitude born at generator order ``n`` with grade ``n/2``. Known products ``K_n^<``, Gram coefficients, and the active flag are computed from it.
- `native_static_solve` is the representation-free core of one slot, in Kossakowski coordinates of a fixed frame. It projects onto the dark quotient, takes the full affine PSD section of ``\Delta_n + \operatorname{Ran}\Phi_n``, factors ``P_n`` into newborn channels, and lifts the bright remainder onto corrections of active amplitudes.
- `native_static_step` runs one slot over a `NativeRepresentation`, which supplies the Kossakowski map, the Hamiltonian read-off, GKSL assembly, and the static gauge directions.
- `native_recurrence` drives the slot through the `BlochFeshbach` and `HoriDeprit` recurrences in the intrinsic total-kick gauge, so both algorithms produce the same intrinsic effective coefficients.
- `DenseLiouvilleRepresentation` is the finite-dimensional backend and the certified reference.

**The canonical section on charged sectors is the resolvent section, not the minimum-norm one.** A `HomologicalInverse` returns a static gauge with ``[\bar{\mathcal L}, S] = -\Pi_{\rm reg}\,\mathcal D[\pi_n\Delta_n\pi_n^\dagger]`` on the sectors where ``\operatorname{ad}_{\bar{\mathcal L}}`` is regular, together with the gauge directions of the remaining sector, where the affine PSD section applies. `NoHomologicalInverse` has no regular sector and reproduces the minimum-norm section everywhere. For a U(1)-symmetric ``\bar{\mathcal L}`` the regular sectors are the charged ones. There ``\Phi_n`` is block diagonal in the superoperator charge, the PSD constraint is inactive, and the resolvent section differs from minimum-norm only by a static gauge in ``\ker\Phi_n``. The two finite generators are therefore similar up to ``O(\varepsilon^{N+1})``. The resolvent section needs no metric on the gauge algebra, so it exists in localized or symbolic coefficient algebras where minimum-norm does not. For driven Kerr it is what the triangular recycling ladder ``(D+J)^{-1}`` computes.

## Public gauge

The construction is exposed as its own gauge, `GKSLNormalForm(; algorithm)`, not as a mode of `VanVleck`. A different gauge is a different struct: `VanVleck` keeps its exact meaning ``\langle\mathcal K\rangle = 0``, while the GKSL normal form adds the static similarity ``W = e^{S}`` and is therefore not zero-average whenever some ``S_n \neq 0``. For Liouvillian generators `floquet_expansion(system, GKSLNormalForm(), order)` returns a `FloquetExpansion` whose completion slot holds a `NativeRealization`, so `effective_generator`, `hamiltonian`, `channels`, `kossakowski`, and `dissipative_frame` behave as after `Gram()`, and `positive_completion` refuses the already GKSL result. For Hamiltonian generators the gauge returns the van Vleck coefficients, since every ``S_n`` vanishes.

## Exact cutoff-free representation

The production backend is `AlgebraicLiouvilleRepresentation` over an exact `OperatorAlgebra` of boson, phase-space, and finite-level sites. `SQALowering` maps Fock, position and momentum, NLevel, and Pauli operators onto it and lifts results back to SQA expressions. Bosonic modes are never truncated. Coefficients are exact rationals; the driver runs in checked `Rational{Int128}` and repeats the expansion in `Rational{BigInt}` only on `OverflowError`. The frame of monomials grows on demand.

Finite Fock truncations are not a substitute. In a truncated space ``[a, a^\dagger] = 1 - d\,|d-1\rangle\langle d-1|``, and for driven Kerr with loss the order-two dark target is then localized on the top Fock levels with a norm that grows with the cutoff. The dense backends remain certified oracles for finite algebras and for the truncated models themselves.

## Positive static section

The dark equation ``\Delta_n + \Phi_n(S_n) = P_n \succeq 0`` is solved in exact arithmetic without a numerical optimizer:

1. Facial reduction: if a principal block of ``P_n`` is untouched by every remaining gauge direction, positivity forces the rows over its kernel to vanish. These linear constraints are imposed until none remain. An indefinite fixed block proves that the gauge family has no PSD lift.
2. Canonical cancellation: within the constrained affine family the least-squares point in the frame metric, with minimum gauge norm, is taken. If it is not PSD, rows of ``P_n`` are cancelled in frame order wherever the constraints stay consistent, and the least-squares point is recomputed.
3. The result is certified by an exact rate-weighted LDL factor, which also returns the newborn channels ``P_n = \sum_j d_j l_j l_j^\dagger`` with rational rates.

The static gauge is polynomial. Per order the smallest gauge degree that admits a certified lift is used, so the section is unique given the frame order. This convention changes only ``O(\varepsilon^{N+1})`` data and never the spectrum through the retained order.

For driven Kerr with static loss the native coefficients equal the van Vleck coefficients through order three, a single channel ``a^2 - 2a^\dagger a`` is born at order four, order five adds none, and order six has no polynomial lift for static gauges up to degree three. That is the polynomial obstruction of the derivation note; resolving it requires the localized coefficient algebra and is not implemented.

## Considered options

- **Van Vleck followed by positive completion.** Rejected as the production architecture: completion repairs a fixed retained generator, while a native expansion must be free to change the retained representative and compensate in the micromotion.
- **Minimum-norm section on every sector.** Rejected for charged sectors because it depends on a Frobenius metric on the gauge algebra that a localized algebra does not carry. It stays the section of the neutral block, where the affine PSD constraint lives.
- **Exact cancellation of the whole charged residual.** Rejected because it gauges away charged Hamiltonian terms and breaks ``H_n^{\rm native} = H_n^{\rm ordinary}`` when there are no jumps.
- **Numerical semidefinite programming for the PSD lift.** Rejected because the result would not be exact; a solver remains a diagnostic only.
- **Truncated Fock spaces for bosonic modes.** Rejected because truncation births cutoff-localized channels.
- **Dense Liouville matrices only in the tests.** Rejected because the generic step needs at least one concrete representation for the interface to be checked, and the dense backend is the certified oracle every symbolic backend is compared against.

## Consequences

The public entry point is `GKSLNormalForm`. The research oracles under `test/gksl/native_*` are thin wrappers around the `src/native` driver. The independent order-2 BF oracle and the HD Lie-transform defect checks stay in the tests as separate implementations.

`native_static_solve` uses floating-point eigen- and Dykstra-based geometry. `native_exact_static_solve` is its exact counterpart for the canonical branch, generic over the scalar field. It takes the operator and gauge metrics explicitly, because an exact frame is not orthonormal. It births channels in rate-weighted form ``P_n = \sum_j d_j l_j l_j^\dagger`` by Hermitian congruence, so no square roots appear, and it rejects an indefinite canonical residual rather than projecting onto the PSD cone. A symbolic representation over `Liouvillian` still needs two further pieces: a finite family of static gauge directions generated from the problem, and a localized coefficient algebra ``\mathfrak A_{\rm poly}[d_q^{-1}]`` for the homological inverse. Physical resonances, zeros of a Bohr multiplier on the physical spectrum, must be rejected rather than localized.

The charged section supersedes the minimum-norm convention for charged sectors in item 3 of Theorem 12 of the native GKLS derivation note.

## Gate

`make test` runs the testsets that hold this decision:

- `test/native/recurrence.jl`: "native recurrence: BF and HD agree through order 4"
- `test/native/static_solve.jl`: "native static slot uses the full affine PSD slice"
- `test/native/exact_static_solve.jl`: "exact static slot matches the float core in orthonormal coordinates" and "exact static slot is covariant under a non-orthonormal frame"
- `test/gksl/native_kerr_fock.jl`: "native GKLS BF/HD finite-Fock driven Kerr"
- `test/gksl/native_kerr_graded_slot.jl`: "graded native slot through order 4 on driven Kerr" and "graded native slot keeps the Hamiltonian branch" (finite-Fock truncated model)
- `test/native/gksl_normal_form.jl`: "GKSLNormalForm coincides with VanVleck wherever no static gauge is needed" and "GKSLNormalForm driven Kerr is GKSL with a drive-induced channel"
- `test/native/sqa_driver.jl`: "cutoff-free driven Kerr births one exact channel at order four"

`make jet` holds optimizer stability of the native layer in `test/quality/JET.jl`.

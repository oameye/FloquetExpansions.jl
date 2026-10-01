# The native GKLS expansion solves its jump amplitudes inside the BF/HD recurrence

The generator-native GKLS expansion does not compute an ordinary Liouvillian Van Vleck expansion and repair its Kossakowski matrix afterwards. At every order the static slot of the Bloch–Feshbach or Hori–Deprit recurrence is chosen so that the order-``n`` coefficient already has the graded Gram form ``c(E_n) = K_n^< + K_n^{\rm lin} + P_n``, and the accepted ``E_n`` is fed back into the same recurrence. The static gauge ``S_n``, the coherent part ``H_n``, and the graded jump amplitudes are therefore state of the recurrence. Every finite truncation is GKLS by construction, and positive completion (ADR 0009) is never called.

The implementation lives in `src/native/`, after `gksl/` and independent of `completion/`:

- `GradedChannel{T}` holds an amplitude born at generator order ``n`` with grade ``n/2``. Known products ``K_n^<``, Gram coefficients, and the active flag are computed from it.
- `native_static_solve` is the representation-free core of one slot, in Kossakowski coordinates of a fixed frame. It projects onto the dark quotient, takes the full affine PSD section of ``\Delta_n + \operatorname{Ran}\Phi_n``, factors ``P_n`` into newborn channels, and lifts the bright remainder onto corrections of active amplitudes.
- `native_static_step` runs one slot over a `NativeRepresentation`, which supplies the Kossakowski map, the Hamiltonian read-off, GKSL assembly, and the static gauge directions.
- `native_recurrence` drives the slot through the `BlochFeshbach` and `HoriDeprit` recurrences in the intrinsic total-kick gauge, so both algorithms produce the same intrinsic effective coefficients.
- `DenseLiouvilleRepresentation` is the finite-dimensional backend and the certified reference.

**The canonical section on charged sectors is the resolvent section, not the minimum-norm one.** A `HomologicalInverse` returns a static gauge with ``[\bar{\mathcal L}, S] = -\Pi_{\rm reg}\,\mathcal D[\pi_n\Delta_n\pi_n^\dagger]`` on the sectors where ``\operatorname{ad}_{\bar{\mathcal L}}`` is regular, together with the gauge directions of the remaining sector, where the affine PSD section applies. `NoHomologicalInverse` has no regular sector and reproduces the minimum-norm section everywhere. For a U(1)-symmetric ``\bar{\mathcal L}`` the regular sectors are the charged ones. There ``\Phi_n`` is block diagonal in the superoperator charge, the PSD constraint is inactive, and the resolvent section differs from minimum-norm only by a static gauge in ``\ker\Phi_n``. The two finite generators are therefore similar up to ``O(\varepsilon^{N+1})``. The resolvent section needs no metric on the gauge algebra, so it exists in localized or symbolic coefficient algebras where minimum-norm does not. For driven Kerr it is what the triangular recycling ladder ``(D+J)^{-1}`` computes.

## Considered options

- **Van Vleck followed by positive completion.** Rejected as the production architecture: completion repairs a fixed retained generator, while a native expansion must be free to change the retained representative and compensate in the micromotion.
- **Minimum-norm section on every sector.** Rejected for charged sectors because it depends on a Frobenius metric on the gauge algebra that a localized algebra does not carry. It stays the section of the neutral block, where the affine PSD constraint lives.
- **Exact cancellation of the whole charged residual.** Rejected because it gauges away charged Hamiltonian terms and breaks ``H_n^{\rm native} = H_n^{\rm ordinary}`` when there are no jumps.
- **Dense Liouville matrices only in the tests.** Rejected because the generic step needs at least one concrete representation for the interface to be checked, and the dense backend is the certified oracle every symbolic backend is compared against.

## Consequences

The native expansion is internal and has no public entry point yet. The research oracles under `test/gksl/native_*` are thin wrappers around the `src/native` driver. The independent order-2 BF oracle and the HD Lie-transform defect checks stay in the tests as separate implementations.

`native_static_solve` uses floating-point eigen- and Dykstra-based geometry. A symbolic representation over `Liouvillian` therefore needs three further pieces before it can run: an exact static solve for the canonical branch (``\Delta_n \in -\operatorname{Ran}\Phi_n``, with Hermitian congruence in place of eigendecomposition), a finite family of static gauge directions generated from the problem, and a localized coefficient algebra ``\mathfrak A_{\rm poly}[d_q^{-1}]`` for the homological inverse. Physical resonances, zeros of a Bohr multiplier on the physical spectrum, must be rejected rather than localized.

The charged section supersedes the minimum-norm convention for charged sectors in item 3 of Theorem 12 of the native GKLS derivation note.

## Gate

`make test` runs the testsets that hold this decision:

- `test/native/recurrence.jl`: "native recurrence: BF and HD agree through order 4"
- `test/native/static_solve.jl`: "native static slot uses the full affine PSD slice"
- `test/gksl/native_kerr_fock.jl`: "native GKLS BF/HD finite-Fock driven Kerr"
- `test/gksl/native_kerr_graded_slot.jl`: "graded native slot through order 4 on driven Kerr" and "graded native slot keeps the Hamiltonian branch"

`make jet` holds optimizer stability of the native layer in `test/quality/JET.jl`.

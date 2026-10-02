# The GKSL normal form solves its jump amplitudes inside the BF/HD recurrence

The generator-native GKSL expansion does not compute an ordinary Liouvillian Van Vleck expansion and repair its Kossakowski matrix afterwards. At every order the static slot of the Bloch–Feshbach or Hori–Deprit recurrence is chosen so that the order-``n`` coefficient already has the graded Gram form ``c(E_n) = K_n^< + K_n^{\rm lin} + P_n``, and the accepted ``E_n`` is fed back into the same recurrence. The static gauge ``S_n``, the coherent part ``H_n``, and the graded jump amplitudes are therefore state of the recurrence. Every finite truncation is GKSL by construction, and positive completion (ADR 0009) is never called.

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

The production backend is `AlgebraicLiouvilleRepresentation` over an exact `OperatorAlgebra` of boson, phase-space, finite-level, and spin sites. A spin site is the size-independent universal enveloping algebra of su(2) in the PBW basis S₊^a S_z^b S₋^c, stored in three monomial slots where every other site uses two. `SQALowering` maps Fock, position and momentum, NLevel, Pauli, and Spin operators onto it and lifts results back to SQA expressions. Bosonic modes are never truncated. Coefficients are exact rationals; the driver runs in checked `Rational{Int128}` and repeats the expansion in `Rational{BigInt}` only on `OverflowError`. The frame of monomials grows on demand.

Finite Fock truncations are not a substitute. In a truncated space ``[a, a^\dagger] = 1 - d\,|d-1\rangle\langle d-1|``, and for driven Kerr with loss the order-two dark target is then localized on the top Fock levels with a norm that grows with the cutoff. The dense backends remain certified oracles for finite algebras and for the truncated models themselves.

## Positive static section

The dark equation ``\Delta_n + \Phi_n(S_n) = P_n \succeq 0`` is solved in exact arithmetic without a numerical optimizer:

1. Facial reduction: if a principal block of ``P_n`` is untouched by every remaining gauge direction, positivity forces the rows over its kernel to vanish. These linear constraints are imposed until none remain. An indefinite fixed block proves that the gauge family has no PSD lift.
2. Canonical cancellation: within the constrained affine family the least-squares point in the frame metric, with minimum gauge norm, is taken. If it is not PSD, rows of ``P_n`` are cancelled in frame order wherever the constraints stay consistent, and the least-squares point is recomputed.
3. The result is certified by an exact rate-weighted LDL factor, which also returns the newborn channels ``P_n = \sum_j d_j l_j l_j^\dagger`` with rational rates.

The static gauge is polynomial. Per order the smallest gauge degree that admits a certified lift is used, so the section is unique given the frame order. This convention changes only ``O(\varepsilon^{N+1})`` data and never the spectrum through the retained order.

For driven Kerr with static loss the native coefficients equal the van Vleck coefficients through order three, a single channel ``a^2 - 2a^\dagger a`` is born at order four, and order five adds none. Order six has no polynomial static gauge: the polynomial families of degree up to five all fail, with an indefinite fixed block or inconsistent kernel rows. Only the last retained order can be rescued, by the terminal virtual gauge below; resolving orders above the first non-polynomial gauge needs the localized coefficient algebra and is not implemented.

### Terminal virtual gauge

Every nonzero difference ``\lambda_\mu - \lambda_\nu`` of the averaged generator has modulus at least ``\min(\kappa/2, 2\Delta)``, so there is no resonance. On the Bargmann lattice ``e_{m,m'} = |m)(m'|`` with ``|m) = a^{\dagger m}|0\rangle``, ``\operatorname{ad}_{\mathcal L_0}`` is a diagonal Bohr multiplier plus the loss coupling ``\kappa\, m m'\, e_{m-1,m'-1}``. Charged sectors lie in its range, but the ladder solution of ``[\mathcal L_0, S] = Y`` never terminates, so the charged ``S_6`` is an infinite sum of monomials, in no finite polynomial algebra and in no finite ring of rational functions of ``n``. Orders four and five are polynomial only because their ladders terminate. The order-four channel ``a^2 - 2a^\dagger a`` mixes charges ``-2`` and ``0``, and positivity at order six ties its charged cross entries to neutral ones, which no polynomial gauge can produce.

At the last retained order ``N`` the gauge never has to be stored. The neutral sector of ``Y`` is solvable exactly when the lattice invariants ``\ell_{m,m'}(Y) = \langle L_{m,m'}|Y|R_{m,m'}\rangle`` vanish, where ``L`` and ``R`` are the left and right eigenvectors of ``\operatorname{ad}_{\mathcal L_0}`` on the lattice. When every polynomial gauge family raises `NativePositivityError` at ``n = N``, the recurrence retries with the virtual family: the neutral polynomial parameters of degree up to ``N`` are restricted to the exact nullspace of the invariants on a lattice window, charged Kossakowski entries are cancelled as images, and the cross units between the charges mixed by the active channels enter as ordinary gauge images of `native_exact_static_solve`. The accepted ``E_N = \hat V_N + Y`` is a polynomial GKSL generator with ``\ell(Y) = 0`` on a window far beyond the fitted one, so its quasi-spectrum agrees with van Vleck through order ``N``. For driven Kerr this births ``a^\dagger`` with rate ``18/78125`` and ``a^{\dagger 2} a^2`` with rate ``5832/78125``, both ``\propto \omega_d^{-6}``, for Bloch–Feshbach and Hori–Deprit alike.

The step carries no stored ``S_N``: `NativeStaticStep.virtual` is set, ``S_N`` is recorded as zero, and the order is listed in `NativeRealization.virtual_orders`. The effective generator, Hamiltonian, channels and Kossakowski matrices are exact. The micromotion of a virtual order keeps its oscillatory harmonics and omits the static one, which is the formal operator ``S_N``.

Limits: the virtual gauge needs one bosonic site, an averaged generator that is number-diagonal plus loss ``a\rho a^\dagger`` without resonances, and the last retained order. A polynomial failure below the last order throws an `ArgumentError` that says an explicit static gauge is required; an unsupported model rethrows the original `NativePositivityError` with the reason appended. The invariants are fitted on a finite window and certified on a larger one, not proved for all lattice points. Models with symbolic parameters treat a virtual order as a nonzero gauge and still require numeric values.

## Symbolic parameters

The exact section runs over numbers. For a model with symbolic parameters `GKSLNormalForm` first runs the exact expansion at two generic rational points, ratios of distinct primes, to decide the structure. If every static gauge vanishes there, the native expansion is the graded positive completion of the van Vleck expansion, so the symbolic result is `positive_completion` of the symbolic van Vleck expansion with `Gram()`, carrying its positivity conditions on the parameters.

Otherwise the exact solver is treated as a black box and its outputs are reconstructed as rational functions of the parameters. The first structure point is the reference. Positivity tests, pivot choices, facial reductions and births are sign decisions, so the outputs are piecewise rational and only the sign region of the reference point is sampled. Every sample must reproduce the reference structure: the monomial frame, channel onsets, virtual and static-gauge orders, and the set of nonzero coefficients. A sample that does not throws an `ArgumentError` instead of mixing two branches. The reconstructed outputs are every scalar the lifting reads: the effective components, the Hamiltonian components, the micromotion harmonics, the channel weights and retained amplitudes, the retained and graded Kossakowski entries, and the finite-generator terms.

The reconstruction runs in three steps over ``\mathbb{Q}(i)``, with small-height sample points so that the solver stays on its `Int128` path. Univariate rational reconstruction with early termination along one line per parameter gives the numerator and denominator degree in each parameter. Along the ray ``t\,x_0`` a polynomial output gives its lowest and highest total degree, which captures the homogeneity that the ``\omega^{-n}`` scaling imposes; an output with a denominator uses a shifted ray for its total degrees instead. A linear solve over the resulting monomial support then fixes the coefficients, adding random sample points until the system has full rank. Before the result is returned it is checked exactly at the reference point and at two fresh points. Each reference value is then mapped to its rational function, and `native_floquet_expansion` lifts the reference data through that map, so the symbolic and numeric paths share the lifting code. Two different functions with the same reference value make the map ambiguous and throw.

The positivity conditions of the reconstructed result are the parameter signs of the sampled region and the channel rates. The result is certified on that region only; other sign regions may carry other rational functions. A single run over the field ``\mathbb{Q}(i)(p_1,\dots,p_m)`` would avoid interpolation and is faster at Kerr order five; issue #389 tracks that route through Nemo.

Exact inputs must reach the expansion exactly. SQA stores a native factor of a raw symbolic coefficient as `ComplexF64`, so the `-i` of a Hamiltonian action would turn the `2//5` of `(2//5) cos(\omega t)` into a float before `harmonics` splits the phase. `Liouvillian` therefore promotes Gaussian-integer native factors to exact constants before they meet a raw coefficient. The native lowering still accepts floats through `rationalize`, but a float input is not an exact model.

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

The charged section supersedes the minimum-norm convention for charged sectors in item 3 of Theorem 12 of the native GKSL derivation note.

## Gate

`make test` runs the testsets that hold this decision:

- `test/native/recurrence.jl`: "native recurrence: BF and HD agree through order 4"
- `test/native/static_solve.jl`: "native static slot uses the full affine PSD slice"
- `test/native/exact_static_solve.jl`: "exact static slot matches the float core in orthonormal coordinates" and "exact static slot is covariant under a non-orthonormal frame"
- `test/gksl/native_kerr_fock.jl`: "native GKSL BF/HD finite-Fock driven Kerr"
- `test/gksl/native_kerr_graded_slot.jl`: "graded native slot through order 4 on driven Kerr" and "graded native slot keeps the Hamiltonian branch" (finite-Fock truncated model)
- `test/native/gksl_normal_form.jl`: "GKSLNormalForm coincides with VanVleck wherever no static gauge is needed" and "GKSLNormalForm driven Kerr is GKSL with a drive-induced channel"
- `test/native/sqa_driver.jl`: "cutoff-free driven Kerr births one exact channel at order four"
- `test/native/gksl_virtual_gauge.jl`: the order-six newborn channels, the lattice invariants of the virtual defect with positive and negative controls, BF/HD agreement, positivity, and the order restriction
- `test/native/gksl_symbolic.jl`: "symbolic GKSLNormalForm without a static gauge is the Gram completion", "symbolic GKSLNormalForm reconstructs the order-four static gauge" and "symbolic GKSLNormalForm with some parameters exact"
- `test/native/rational_reconstruction.jl`: "exact rational reconstruction recovers sparse rational functions" and "exact rational reconstruction rejects a point that breaks the structure"

`make jet` holds optimizer stability of the native layer in `test/quality/JET.jl`.

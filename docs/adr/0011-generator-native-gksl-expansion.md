# The GKSL normal form solves its jump amplitudes inside the BF/HD recurrence

The generator-native GKSL expansion does not compute an ordinary Liouvillian Van Vleck expansion and repair its Kossakowski matrix afterwards. At every order the static slot of the Bloch–Feshbach or Hori–Deprit recurrence is chosen so that the order-``n`` coefficient already has the graded Gram form ``c(E_n) = K_n^< + K_n^{\rm lin} + P_n``, and the accepted ``E_n`` is fed back into the same recurrence. The static gauge ``S_n``, the coherent part ``H_n``, and the graded jump amplitudes are therefore state of the recurrence. Every finite truncation is GKSL by construction, and positive completion (ADR 0009) is never called.

The implementation lives in `src/native/`, after `gksl/` and independent of `completion/`:

- `GradedChannel{T}` holds an amplitude born at generator order ``n`` with grade ``n/2``. Known products ``K_n^<``, Gram coefficients, and the active flag are computed from it.
- `native_exact_static_solve` is the representation-free core of one slot, in Kossakowski coordinates of a fixed frame and exact in the scalar field. It projects onto the dark quotient, takes an affine positive-semidefinite section of ``\Delta_n + \operatorname{Ran}\Phi_n``, factors ``P_n`` into newborn channels, and lifts the bright remainder onto corrections of active amplitudes.
- `native_static_step` runs one slot over a `NativeRepresentation`, which supplies the Kossakowski map, the Hamiltonian read-off, GKSL assembly, and the static gauge directions.
- `native_recurrence` drives the slot through the `BlochFeshbach` and `HoriDeprit` recurrences in the intrinsic total-kick gauge, so both algorithms produce the same intrinsic effective coefficients.
- `AlgebraicLiouvilleRepresentation` is the only `NativeRepresentation`. It works over an exact operator algebra, and `SQALowering` maps SQA operators onto it and lifts results back.

## Public gauge

The construction is exposed as its own gauge, `GKSLNormalForm(; algorithm)`, not as a mode of `VanVleck`. A different gauge is a different struct: `VanVleck` keeps its exact meaning ``\langle\mathcal K\rangle = 0``, while the GKSL normal form adds the static similarity ``W = e^{S}`` and is therefore not zero-average whenever some ``S_n \neq 0``. For Liouvillian generators `floquet_expansion(system, GKSLNormalForm(), order)` returns a `FloquetExpansion` whose completion slot holds a `NativeRealization`, so `effective_generator`, `hamiltonian`, `channels`, `kossakowski`, and `dissipative_frame` behave as after `Gram()`, and `positive_completion` refuses the already GKSL result. For Hamiltonian generators the gauge returns the van Vleck coefficients without a completion state, since every ``S_n`` vanishes.

## Exact cutoff-free representation

The operator algebra covers boson, phase-space, finite-level, and spin sites. A spin site is the size-independent universal enveloping algebra of su(2) in the PBW basis S₊^a S_z^b S₋^c, stored in three monomial slots where every other site uses two. `SQALowering` maps Fock, position and momentum, NLevel, Pauli, and Spin operators onto it. Bosonic modes are never truncated. The lowering accepts only exact `Integer` and `Rational` coefficients with no free symbol, and throws an `ArgumentError` for a floating-point value, for a symbolic parameter, and for two operators of different names on one Hilbert space. Symbolic parameters reach it only as exact values, through the reconstruction described below. The driver runs in checked `Rational{Int128}` and repeats the expansion in `Rational{BigInt}` only on `OverflowError`. The frame of monomials grows on demand. An averaged generator that is not of GKSL form, that is not trace preserving, throws an `ArgumentError`.

Finite Fock truncations are not a substitute. In a truncated space ``[a, a^\dagger] = 1 - d\,|d-1\rangle\langle d-1|``, and for driven Kerr with loss the order-two dark target is then localized on the top Fock levels with a norm that grows with the cutoff.

Exact inputs must reach the expansion exactly. SQA 0.13 keeps coefficients exact by type (ADR 0007), so the `2//5` of `(2//5) cos(\omega t)` stays rational through the `-i` of a Hamiltonian action and `harmonics` splits the phase without a float intermediate.

## Positive static section

The dark equation ``\Delta_n + \Phi_n(S_n) = P_n \succeq 0`` is solved in exact arithmetic without a numerical optimizer:

1. Facial reduction: if a principal block of ``P_n`` is untouched by every remaining gauge direction, positivity forces the rows over its kernel to vanish. These linear constraints are imposed until none remain. An indefinite fixed block proves that the gauge family has no PSD lift.
2. Canonical cancellation: within the constrained affine family the least-squares point in the frame metric, with minimum gauge norm, is taken. If it is not PSD, rows of ``P_n`` are cancelled in frame order wherever the constraints stay consistent, and the least-squares point is recomputed.
3. The result is certified by an exact rate-weighted LDL factor, which also returns the newborn channels ``P_n = \sum_j d_j l_j l_j^\dagger`` with rational rates, by Hermitian congruence and so without square roots.

The static gauge is polynomial. The gauge directions are tried by increasing degree, and the first family that admits a certified lift is used, so the section is unique given the frame order. This convention changes only ``O(\varepsilon^{N+1})`` data and never the spectrum through the retained order. A slot that no family can lift raises `NativePositivityError`, which `floquet_expansion` reports as an order that needs a static gauge that is not polynomial.

For driven Kerr with static loss the native coefficients equal the van Vleck coefficients through order three, a single channel ``a^2 - 2a^\dagger a`` is born at order four, and order five adds none. Order six has no polynomial static gauge, so the expansion is available through retained order five and `floquet_expansion(..., 7)` throws.

## Considered options

- **Van Vleck followed by positive completion.** Rejected as the production architecture: completion repairs a fixed retained generator, while a native expansion must be free to change the retained representative and compensate in the micromotion.
- **Exact cancellation of the whole charged residual.** Rejected because it gauges away charged Hamiltonian terms and breaks ``H_n^{\rm native} = H_n^{\rm ordinary}`` when there are no jumps.
- **Numerical semidefinite programming for the PSD lift.** Rejected because the result would not be exact; a solver remains a diagnostic only.
- **Truncated Fock spaces for bosonic modes.** Rejected because truncation births cutoff-localized channels.

## Consequences

The public entry point is `GKSLNormalForm`. The independent order-2 BF oracle, the finite-dimensional Floquet-exponent check, and the Gram-completion comparison stay in the tests as separate implementations.

`native_exact_static_solve` takes the operator and gauge metrics explicitly, because an exact frame is not orthonormal, and it rejects an indefinite canonical residual rather than projecting onto the PSD cone.

## Symbolic parameters

A generator with free symbols is expanded by reconstruction, not by symbolic arithmetic in the solver. The solver stays numeric and exact. `symbolic_parameters.jl` evaluates it at rational points of a declared sign region, records every exact output that the lift to SQA reads, and `polynomial_reconstruction.jl` rebuilds each output as a sparse polynomial in the parameters. Outputs that agree at the first sample points share one fit, and each class is fitted at its own support, found from degrees along coordinate lines and a ray. The linear algebra is exact over `Rational{BigInt}`, with real and imaginary parts as separate outputs. Only polynomials are reconstructed: an output that is rational in the parameters fails with an `ArgumentError`, because no model tried so far needed one and a rational fit would add an unexercised route.

The region is part of the result. Every parameter is positive except those the gauge lists in `negative`, and the drive frequency is positive, because an odd onset carries a factor ``\omega^{-n}`` that changes the sign of a rate. The reference point is a fixed ratio of primes in the region, and the exact outputs there define the structure: the frame, the channel onsets, the orders carrying a static gauge, and the set of nonzero outputs.

The reconstruction is certified by sampling, not by proof. After the fit, every output is compared exactly with the solver at `verification` fresh points, each parameter scaled independently over a factor of 729 around values between ``1/9`` and ``9/2``, and each point must reproduce the structure of the reference point. This covers exactly the claim in the docstring: the returned expression equals the exact numeric expansion at the listed points, which `factorization` returns. It does not cover other points of the region. A branch of the solver confined to an unsampled part of the region is invisible, and a reconstruction whose structure changes at a sampled point fails with an `ArgumentError` instead of returning a formula for one piece. Fractional powers of a parameter are rejected when the parameters are collected, since an irrational root has no exact value to sample, and at most six parameters are accepted because the sampling cost grows with their number.

The positivity conditions are the signed parameters, the drive frequency, and every channel rate with its power of the drive frequency. The exact route over a fraction field, which would replace sampling by a proof, is tracked in issue 389.

## Deferred

A terminal virtual gauge for an order whose static gauge is not polynomial is left to follow-up work. The expansion refuses such input rather than approximating it.

## Gate

`make test` runs the testsets that hold this decision:

- `test/native/gksl_normal_form.jl`: "GKSLNormalForm coincides with VanVleck wherever no static gauge is needed", "GKSLNormalForm driven Kerr is GKSL with a drive-induced channel" (which also compares BF with HD), "GKSLNormalForm on a driven oscillator in position and momentum" and "GKSLNormalForm rejects what it cannot expand exactly"
- `test/native/gksl_gram_regression.jl`: "driven Kerr with loss: GKSLNormalForm equals Gram completion of Van Vleck" and "driven qubit with decay: GKSLNormalForm equals Gram completion of Van Vleck"
- `test/native/gksl_finite_oracle.jl`: "GKSLNormalForm spectrum matches the exact Floquet exponents" and "GKSLNormalForm Kossakowski matrix is positive semidefinite at finite frequency"
- `test/native/gksl_phase_independence.jl`: "GKSLNormalForm is independent of the drive phase"
- `test/native/exact_static_solve.jl`: "exact static slot solves the orthonormal fixture", "exact static slot is covariant under a non-orthonormal frame", "exact static slot is covariant under a gauge reparameterization" and "exact static slot rejects an indefinite canonical residual"
- `test/native/graded_channels.jl`: "graded channel Gram coefficients are exact Cauchy products" and "graded channels activate after birth and stay prefix complete"
- `test/native/sqa_driver.jl`: "native expansion equals Van Vleck wherever no static gauge is needed", "cutoff-free driven Kerr births one exact channel at order four", "SQA lowering is exact across Fock, level and Pauli sites" and "the exact driver falls back to BigInt scalars when Int128 overflows"
- `test/native/gksl_symbolic.jl`: the symbolic result equals the numeric expansion at its verification points, equals the Gram completion of the van Vleck expansion where no static gauge is needed, reproduces the order-four channel of driven Kerr, reports its region, and rejects what it cannot certify
- `test/native/polynomial_reconstruction.jl`: the interpolation recovers sparse polynomials per sign region, is bound to the declared region, and rejects a branch beyond the fitted magnitudes and a coefficient that is not a polynomial
- `test/native/spin_sites.jl`: "native expansion of a driven spin equals Van Vleck" and "GKSLNormalForm on a driven spin is completely positive"

`make jet` holds optimizer stability of the graded channels in `test/quality/JET.jl`, "native graded channel optimizer stability". The exact static solve is not covered.

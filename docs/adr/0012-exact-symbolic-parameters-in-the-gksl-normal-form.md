# The GKSL normal form expands symbolic parameters exactly over a rational-function field

ADR 0011 left symbolic parameters to follow-up work, and the native lowering refused them. A model is almost always written with parameters, so `GKSLNormalForm` now expands a generator whose coefficients are polynomials in symbolic parameters with exact rational coefficients. The native recurrence runs unchanged over the scalar field ``\mathbb Q(p_1, \dots, p_k)``. Its rational-function elements are Nemo's `FracFieldElem{QQMPolyRingElem}`, backed by FLINT, so every zero test is exact and every element has one canonical form.

The implementation lives in `src/native/parameter_field.jl`:

- `ParameterScalar <: Real` wraps an element of the fraction field of the cached Nemo polynomial ring in at most eight variables. The native solver uses `Complex{ParameterScalar}` as its scalar, as it uses `Complex{Rational{Int128}}` for numeric input.
- `native_coefficients` reads the parameters of a generator and returns `ExactCoefficients()` when there are none, so numeric input takes the existing checked-`Int128` path without change.
- A parameter declared negative is lowered as ``p = -x``. Every field variable is then positive.
- SQA splits a complex parameter ``\alpha`` into the real symbols `real(α)` and `imag(α)`, which become two field variables.
- The drive frequency may not appear in a coefficient. The recurrence orders terms by explicit powers of ``1/\omega``, so an ``\omega`` inside an order-zero coefficient would silently break that ordering.

## Sign decisions

The exact static solve decides signs only in its rate-weighted ``LDL^\dagger`` certificate. A sign is decided when the numerator and the denominator of the rational function each have coefficients of one sign. Every field variable is positive, so the sign then holds throughout the sign region. A gauge family whose certificate meets any other sign is treated as having no lift throughout the region, and the next family is tried, so the symbolic section is the first family certified on the whole region. Where no family is, the expansion throws an `ArgumentError` that names the expression, for instance `κ - γ`, rather than deciding it at a sample point. For driven Kerr, the driven qubit, and their variants with dephasing or pumping every sign decision is a monomial and no family is skipped.

Pivots are chosen for generic parameters, but no pivot is recorded. On the dense open set where every pivot is nonzero the returned generator is of GKSL form, its rates are nonnegative, and its spectrum agrees with van Vleck. These are closed conditions on the returned coefficients, so they extend by continuity to every point where the coefficients are finite. `regularity_conditions` therefore returns the irreducible factors of the coefficient denominators, for instance ``\kappa`` for driven Kerr with loss and dephasing, and `positivity_conditions` returns the signed parameters and the drive frequency.

## Considered options

- **Sampling at numeric points and rational reconstruction** (draft PR #396). Rejected: reconstruction needs degree bounds that the recurrence does not supply, and a wrong guess returns a plausible but wrong rate.
- **Symbolics `Num` as the scalar.** Rejected: it has no reliable zero test, and the exact static solve branches on zeros.
- **DynamicPolynomials with a hand-written fraction.** Rejected: slower than Nemo on every measured solve, and gcd normalization would be ours to maintain.
- **Nemo as a package extension.** Rejected: symbolic parameters occur in almost every model, so the field is core functionality, not an optional backend.
- **Deciding a sign at a fixed reference point.** Rejected: a sign that changes inside the declared region would silently restrict the result to part of it.

## Consequences

Nemo is a direct dependency, which adds the FLINT binaries to every installation. Its compat follows the range that Symbolics accepts. Driven Kerr with loss expanded to order five in four parameters takes well under a second, against about a quarter of a second numerically.

## Gate

`make test` runs the testsets that hold this decision:

- `test/native/gksl_symbolic.jl`: "symbolic GKSLNormalForm equals the numeric expansion at substituted points", "symbolic GKSLNormalForm states where it holds" and "symbolic GKSLNormalForm rejects what it cannot decide"

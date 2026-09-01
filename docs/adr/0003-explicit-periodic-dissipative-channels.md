# Explicit periodic dissipative channel forms

The `Liouvillian(H; collapse_operators=(), jumps=(), rates=())` constructor accepts complete collapse operators `L`, contributing `D[L]`, and rate-weighted jump channels `(γ, J)`, contributing `γD[J]`. Rates are ordinary symbolic scalar coefficients in the initial implementation.

Periodic dissipation is represented by lowering each Fourier component into a Liouvillian and placing those components in `PeriodicGenerator`; there is no separate time-dependent-rate wrapper. This supports periodic dependence in the complete collapse operator, the scalar rate, or both, while keeping the expansion input uniform.

The package does not add `NonnegativeRate` types, symbolic sign conditions, or matrix-valued rate data.

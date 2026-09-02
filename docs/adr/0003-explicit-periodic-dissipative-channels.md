# Explicit periodic dissipative channel forms

The `Liouvillian(H; channels=())` constructor accepts explicit channel values. `collapse(L)`
contributes `D[L]`, while `jump(J; rate=γ)` contributes `γD[J]`. Rates are ordinary symbolic
scalar coefficients in the initial implementation.

For a time-dependent input, construct the native map first and call
`harmonics(L, w, t)`. The collector decomposes every left and right operator factor and every
scalar coefficient, then places the result in the common `PeriodicGenerator{Liouvillian}`
container. This supports periodic dependence in the Hamiltonian, complete collapse operator,
jump operator, rate, or any combination without introducing a separate periodic-channel type.
The `channels` keyword of `floquet_expansion(H, w, t, gauge, order)` provides the same lowering as
a convenience. The two channel forms are additive for independent processes; they are not inferred
to be the same channel or deduplicated.

The package does not add `NonnegativeRate` types, symbolic sign conditions, or matrix-valued rate data.

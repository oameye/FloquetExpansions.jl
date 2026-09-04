# Explicit dissipative channel forms

`liouvillian(H; channels=())` accepts `collapse(L)` for a complete collapse operator and `jump(J, γ)` for a separate scalar rate, contributing `D[L]` and `γD[J]`. These forms are additive and are not inferred to represent the same channel.

`jump(J, γ)` is a physical rate-weighted channel: `γ` must be real and is asserted nonnegative. Provably negative numeric rates and non-real rates are rejected. A symbolic rate expression is asserted nonnegative as a whole; the package does not infer separate sign assumptions for its constituent symbols. `collapse(L)` remains the representation for a complete collapse operator with amplitudes and phases already folded into `L`.

For time-dependent inputs, `harmonics` decomposes a `Liouvillian`'s operator factors and scalar coefficients into `PeriodicGenerator{Liouvillian}`. The `channels` keyword on `floquet_expansion` performs this conversion automatically, supporting periodic dependence in the Hamiltonian, collapse operator, jump operator, rate, or any combination.

The high-level `floquet_expansion(H, ωd, t, gauge, order; channels=...)` path may retain internal microscopic dissipative provenance for later positive-completion work. Generic `Liouvillian` and `PeriodicGenerator` values remain algebraic and do not carry channel provenance through arbitrary addition, composition, or commutators.

# Explicit dissipative channel forms

`Liouvillian(H; channels=())` accepts `collapse(L)` for a complete collapse operator and `jump(J, γ)` for a separate scalar rate, contributing `D[L]` and `γD[J]`. These forms are additive and are not inferred to represent the same channel.

For time-dependent inputs, `harmonics` decomposes a `Liouvillian`'s operator factors and scalar coefficients into `PeriodicGenerator{Liouvillian}`. The `channels` keyword on `floquet_expansion` performs this conversion automatically, supporting periodic dependence in the Hamiltonian, collapse operator, jump operator, rate, or any combination.

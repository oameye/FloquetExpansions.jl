# Unified Floquet expansion model

Hamiltonian and Liouvillian dynamics use one `FloquetExpansion` result and one `floquet_expansion` function. Its order-resolved state is represented by `generator`, `kick_components`, `kick_derivative_components`, `dressed_generator`, `dressed_kick_derivative`, `effective_components`, `gauge`, and `order`, with the component type shared throughout the homological recursion.

The generic accessors are `effective_generator` and `micromotion`. `effective_generator(F)` reattaches and sums the computed inverse-frequency orders, while `effective_generator(F, n)` returns one order contribution. `micromotion(F)` returns the periodic kick series and `micromotion(F, n)` returns one positive-order kick contribution. The existing Hamiltonian-only accessors remain compatibility aliases, while a parallel `FloquetLiouvillianExpansion` type is not retained.

# One result model and generic accessors

Hamiltonian and Liouvillian dynamics share `FloquetExpansion` and `floquet_expansion`. Consumers use `effective_generator` and `micromotion`; order-specific accessors return the corresponding inverse-frequency contributions, and no separate Liouvillian result type is introduced.

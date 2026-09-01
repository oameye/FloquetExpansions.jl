# SQA-backed left/right Liouvillian representation

A Liouvillian is represented as a sparse sum of elementary maps `ρ ↦ AρB`, with `A` and `B` stored as SQA operator expressions and scalar coefficients collected eagerly in a `Dict{Tuple{QAdd,QAdd},CNum}`. SQA therefore owns operator multiplication, normal ordering, adjoints, and symbolic coefficient algebra; FloquetExpansions only combines equal left/right actions and distributes composition. Cheap algebraic normalization is eager, while expensive whole-expression simplification remains explicit.

`compose(A, B)` means that `B` acts first. For elementary terms it maps `(Aₗ, Aᵣ) ∘ (Bₗ, Bᵣ)` to `(AₗBₗ, BᵣAᵣ)` and multiplies the scalar coefficients; the resulting equal actions are collected immediately.

The coherent and dissipative constructor inputs remain distinct for clarity, but the stored algebra is uniform. This is necessary because composition and commutators produce cross terms that are not generally expressible as separate Hamiltonian and standard dissipator fields. The representation has no custom Liouvillian expression tree or symbolic positivity conditions. Hilbert-space compatibility and product-space structure are inherited from the supplied SQA expressions; FloquetExpansions does not perform a second space-validation layer.

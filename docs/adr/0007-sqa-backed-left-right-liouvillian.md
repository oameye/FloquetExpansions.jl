# SQA-backed left/right Liouvillian representation

Liouvillians use a sparse collected sum of elementary maps `ρ ↦ AρB`, with SQA owning operator multiplication, adjoints, normal ordering, and coefficient algebra. This uniform representation keeps composition and commutators closed even when terms no longer have a Hamiltonian-plus-dissipator presentation; equal actions are collected eagerly, while whole-expression simplification remains explicit.

Term insertion is owned by the Liouvillian module, and Fourier lowering for Liouvillians (`harmonics(::Liouvillian)`) is owned there as well: it consumes module-owned term and harmonic iterators and the `PeriodicGenerator` Fourier seam. The periodic-generator module owns the shared phase → harmonic normalization (`harmonic_index` and its internal phase visitor), while Liouvillian lowering retains its local left/right convolution and sparse term insertion. This keeps sparse storage out of the shared Fourier seam while preserving the eager `Dict` implementation and its measured allocation profile.

The public `terms(L)` iterator exposes the semantic `(left, right, coefficient)` triples without making callers depend on the sparse storage layout. Numerical vectorization remains a separate consumer of this interface rather than part of the symbolic Liouvillian core.

`compose(A, B)` means that `B` acts first and `A` acts second. For elementary actions it maps `(Aₗ, Aᵣ) ∘ (Bₗ, Bᵣ)` to `(AₗBₗ, BᵣAᵣ)`.

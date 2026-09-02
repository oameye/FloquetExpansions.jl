# SQA-backed left/right Liouvillian representation

Liouvillians use a sparse collected sum of elementary maps `ρ ↦ AρB`, with SQA owning operator multiplication, adjoints, normal ordering, and coefficient algebra. This uniform representation keeps composition and commutators closed even when terms no longer have a Hamiltonian-plus-dissipator presentation; equal actions are collected eagerly, while whole-expression simplification remains explicit.

Term insertion is owned by the Liouvillian module, and Fourier lowering consumes module-owned term and harmonic iterators. This keeps sparse storage out of the collector while preserving the eager `Dict` implementation and its measured allocation profile.

`compose(A, B)` means that `B` acts first and `A` acts second. For elementary actions it maps `(Aₗ, Aᵣ) ∘ (Bₗ, Bᵣ)` to `(AₗBₗ, BᵣAᵣ)`.

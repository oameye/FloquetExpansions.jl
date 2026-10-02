# SQA-backed left/right Liouvillian representation

Liouvillians use a sparse collected sum of elementary maps `ρ ↦ AρB`, with SQA owning operator multiplication, adjoints, normal ordering, and coefficient algebra. This uniform representation keeps composition and commutators closed even when terms no longer have a Hamiltonian-plus-dissipator presentation; equal actions are collected eagerly, while whole-expression simplification remains explicit.

Term insertion is owned by the Liouvillian module, and Fourier lowering for Liouvillians (`harmonics(::Liouvillian)`) is owned there as well: it consumes module-owned term and harmonic iterators and the `PeriodicGenerator` Fourier seam. The periodic-generator module owns the shared phase → harmonic normalization (`harmonic_index` and its internal phase visitor), while Liouvillian lowering retains its local left/right convolution and sparse term insertion. This keeps sparse storage out of the shared Fourier seam while preserving the eager `Dict` implementation and its measured allocation profile.

The public `terms(L)` iterator exposes the semantic `(left, right, coefficient)` triples without making callers depend on the sparse storage layout. Numerical vectorization remains a separate consumer of this interface rather than part of the symbolic Liouvillian core.

`compose(A, B)` means that `B` acts first and `A` acts second. For elementary actions it maps `(Aₗ, Aᵣ) ∘ (Bₗ, Bᵣ)` to `(AₗBₗ, BᵣAᵣ)`.

## Amendment: canonical keys

"Equal actions are collected eagerly" originally meant equal `(A, B)` factor pairs. Two representations of one map, such as `(a + a', 1)` and `(a, 1) + (a', 1)`, then stayed apart, and a map equal to zero could hold terms. `iszero` and `==` thus compared representations, not maps: the commutator of `hamiltonian_action(a' * a)` and `dissipator(a)` kept four terms. Term insertion now splits each factor into unit-coefficient monomials after the completeness relation of every `NLevelSpace` is applied, so a key is a pair of basis operators and `iszero` and `==` are map equality up to coefficient canonicalization. Factors with bound symbolic sums keep the whole-factor key, because SQA does not reduce them to a basis. Keys stay `Tuple{QAdd,QAdd}`, so `terms(L)` keeps its triples, now one per monomial pair.


## Amendment: exact coefficients

Splitting factors into monomials moves every factor coefficient into the term coefficient, so the `-i` of a Hamiltonian action now multiplies raw trigonometric coefficients such as `(2//5) cos(ωt)` directly. SQA embeds a native factor of a raw expression as `ComplexF64`, which made such coefficients floats and the van Vleck recurrence of a numeric model inexact. Term insertion, scaling, and composition promote Gaussian-integer native factors to exact constants before they meet a raw coefficient, so exact amplitudes stay exact through `harmonics`.

## Gate

`make test` runs the testsets that hold this decision:

- `test/generators/liouvillian.jl`: "Liouvillian terms expose semantic triples"
- `test/generators/liouvillian.jl`: "Liouvillian arithmetic collects equal terms"
- `test/generators/liouvillian.jl`: "Liouvillian composition is map composition"
- `test/generators/liouvillian.jl`: "zero operator factors produce zero maps"
- `test/generators/liouvillian.jl`: "Liouvillian equality is map equality"
- `test/generators/liouvillian.jl`: "exact drive amplitudes stay exact through harmonics"

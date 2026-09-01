# One periodic-generator expansion core

The van Vleck expansion operates on one `PeriodicGenerator{T}` Fourier container. The prototype stores `components::Dict{Int,T}`, the symbolic drive frequency `wd`, and a concrete `zero_component::T`; zero-valued harmonics are omitted, while the zero component preserves the element type for empty generators and missing-harmonic lookup. `T` may be an SQA operator expression or a symbolic Liouvillian, so Hamiltonian and dissipative systems share one arbitrary-order recursion and one result model.

```@meta
CurrentModule = FloquetExpansions
CollapsedDocStrings = true
```

# GKSL coordinates

A [`Liouvillian`](@ref) is an algebraic left/right superoperator and does not carry a preferred
Lindblad representation. When a dissipative operator frame is known, FloquetExpansions can extract
its coherent Hamiltonian and Kossakowski matrix exactly in the symbolic operator algebra.

A [`DissipativeFrame`](@ref) is an ordered set of linearly independent operator directions modulo
the identity. The ordering is part of the representation because it fixes the coordinates of the
Kossakowski matrix. The frame need not be orthogonal or normalized.

```@docs
DissipativeFrame
```

For a Liouvillian written in a supplied frame ``\{F_\mu\}``, [`kossakowski`](@ref) returns the
Hermitian coefficient matrix ``d`` in

```math
\mathcal L(\rho) = -i[H,\rho]
+ \sum_{\mu\nu} d_{\mu\nu}
\left(F_\mu \rho F_\nu^\dagger
- \frac12\{F_\nu^\dagger F_\mu,\rho\}\right).
```

Extraction uses the two-sided sandwich terms directly and verifies the remaining one-sided action
as a Hamiltonian commutator. No Hilbert-space matrix representation or finite-dimensional cutoff is
introduced.

```@docs
kossakowski
kossakowski_component
hamiltonian
hamiltonian_component
```

For a raw Liouvillian Floquet expansion, Kossakowski coordinates require an explicit frame:

```julia
frame = DissipativeFrame(a, a^2)
d = kossakowski(vv, frame)
d1 = kossakowski_component(vv, frame, 1)
```

The component accessor includes the corresponding inverse-drive-frequency scaling. The coherent
Hamiltonian is defined modulo an additive multiple of the identity, which does not change its
commutator action.

Positive completion and automatic finalized-frame storage are separate opt-in features described by
the CP-preserving completion architecture; raw Floquet expansions remain algebraic by default.

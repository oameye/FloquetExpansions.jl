```@meta
CurrentModule = FloquetExpansions
CollapsedDocStrings = true
```

# System

`FloquetExpansions.jl` describes a driven quantum system through its generator. A closed system is
specified by a Hamiltonian whose commutator gives the generator, while an open system is described
directly by a Liouvillian. The first step for either case is to express the time dependence in a
common Fourier representation.

## Symbolic substrate

FloquetExpansions operates on symbolic operators from
[SecondQuantizedAlgebra.jl](@extref SecondQuantizedAlgebra :doc:`index`). SQA provides Hilbert
spaces, quantum operators, symbolic coefficients, and their algebra; this manual focuses on the
generator and Floquet layers built on top of that substrate.

## Generators of quantum dynamics

For a density operator ``\rho(t)``, a generator ``\mathcal{G}(t)`` defines the equation of motion

```math
\frac{d\rho}{dt} = \mathcal{G}(t)[\rho].
```

### Closed systems: Hamiltonian generators

For a closed system, a Hamiltonian ``H(t)`` induces the commutator action

```math
\mathcal{G}_H(t)[\rho] = -i[H(t),\rho].
```

Use [`hamiltonian_action`](@ref) to turn the Hamiltonian expression into the corresponding
Liouvillian map.

```@docs
hamiltonian_action
```

### Open systems: Liouvillian generators

For an open system, the generator acts directly on density operators:

```math
\frac{d\rho}{dt} = \mathcal{L}(t)[\rho].
```

The [`Liouvillian`](@ref) type is a symbolic linear map represented as a sum of left/right
operator actions. This is an algebraic representation: a general Liouvillian need not be in
Lindblad form and does not, by itself, certify positivity or complete positivity.

```@docs
Liouvillian
```

### GKLS (Lindblad) form

A [Lindblad generator](https://en.wikipedia.org/wiki/Lindbladian), also called a
Gorini–Kossakowski–Sudarshan–Lindblad generator (GKSL or GKLS), is the physically constrained
subclass with the form

```math
\mathcal{L}[\rho]
= -i[H,\rho] + \sum_a \gamma_a\,\mathcal{D}[J_a](\rho),
\qquad \gamma_a \ge 0.
```

The elementary dissipator is

```math
\mathcal{D}[J](\rho)
= J\rho J^\dagger
- \frac{1}{2}\left(J^\dagger J\rho + \rho J^\dagger J\right).
```

The following constructors encode these channel terms. Their docstrings define the distinction
between a complete collapse operator and a bare jump operator with a separate rate.

```@docs
dissipator
```

```@docs
collapse
```

```@docs
jump
```

The package accepts symbolic rates as algebraic coefficients. It does not infer that a symbolic
rate is nonnegative, nor does constructing a [`Liouvillian`](@ref) certify that the result is a
GKLS generator.

### Actions and composition

The left/right representation writes a map as

```math
\mathcal{L}[\rho] = \sum_j c_j A_{j,\mathrm{L}}\rho A_{j,\mathrm{R}}.
```

[`actions`](@ref) exposes the terms as semantic triples `(left, right, coefficient)`. This gives
callers a stable way to inspect a map without depending on its sparse storage.

```@docs
actions
```

Liouvillian maps compose as ordinary linear maps: in `compose(A, B)`, `B` acts first and `A` acts
second. Composition and commutators therefore remain in the same symbolic map algebra.

```@docs
compose
```

## Periodic generators

When the generator is periodic with angular frequency ``\omega_d`` and period
``T_d = 2\pi/\omega_d``, write it as

```math
G(t) = \sum_l G_l\,e^{-i l \omega_d t}.
```

Here `l` is an integer harmonic label and ``G_l`` is an operator or Liouvillian component. A
[`PeriodicGenerator`](@ref) stores these components together with the drive frequency, and its
component type determines the algebra used by the Floquet expansion. Missing harmonics are zero;
the full time-dependent expression is recovered by evaluating the generator at a symbolic time.

```@docs
PeriodicGenerator
```

### Fourier decomposition

Use [`harmonics`](@ref) to decompose a symbolic Hamiltonian into the convention above. The API
also decomposes a [`Liouvillian`](@ref) by applying Fourier decomposition to its operator factors
and scalar coefficients, so time dependence may occur in the Hamiltonian, channels, rates, or any
combination.

```@docs
harmonics
```

### Harmonic support and Fourier calculus

[`support`](@ref) reports the smallest range containing the nonzero harmonic labels.

```@docs
support
```

The Fourier representation turns the basic calculus into operations on harmonic components. The
zeroth component is the period average, and differentiation acts diagonally:

```math
(\partial_t G)_l = -i l\,G_l.
```

```@docs
time_average
```

```@docs
derivative
```

[`antiderivative`](@ref) inverts this operation on a zero-average generator. Its free integration
constant is fixed by the selected [`Gauge`](@ref), and is used by the Floquet expansion.

```@docs
antiderivative
```

Periodic commutators use harmonic convolution, while symbolic simplification is applied component
by component.

```@docs
SecondQuantizedAlgebra.commutator(::PeriodicGenerator{T}, ::PeriodicGenerator{T}) where {T}
```

```@docs
SecondQuantizedAlgebra.simplify(::PeriodicGenerator{T}) where {T}
```

For Hamiltonian components, the adjoint reverses the harmonic label. Hermiticity is consequently
the condition ``G_{-l} = G_l^\dagger`` for every `l`.

```@docs
Base.adjoint(::PeriodicGenerator{SecondQuantizedAlgebra.QAdd})
LinearAlgebra.ishermitian(::PeriodicGenerator{SecondQuantizedAlgebra.QAdd})
```

## Quasienergy representation

[`QuasienergyOperator`](@ref) is an alternative representation in truncated Sambe space. It
assembles the Fourier components into blocks indexed by harmonic labels, which is useful for
quasienergy calculations and for the Floquet–Liouville analogue in open systems. The API docstring
below gives the Hamiltonian and Liouvillian block conventions, truncation, and interpretation of
complex dissipative quasienergies.

```@docs
QuasienergyOperator
```

```@docs
harmonic_range
```

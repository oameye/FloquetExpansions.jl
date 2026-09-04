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

FloquetExpansions does not define a separate operator algebra. It operates on symbolic operators
from [SecondQuantizedAlgebra.jl](@extref SecondQuantizedAlgebra :doc:`index`), or SQA, and adds the
generator and Floquet layers on top of them.

There are two distinct symbolic layers:

- SQA defines Hilbert spaces, quantum operators, and the algebra used when operators are added or
  multiplied.
- [Symbolics.jl](https://symbolics.juliasymbolics.org/stable/) defines scalar symbolic expressions
  such as frequencies, coupling strengths, amplitudes, and time. Use `@variables` to create them;
  declare real quantities with `::Real` and potentially complex quantities with `::Number`.

For complex scalar symbols, prefer `::Number` to `::Complex`.[^complex_symbol]

### Choosing an operator algebra

The Hilbert-space type selects the operator model and its defining relations:

| Physical subsystem | SQA space | Operators | Defining relation |
| --- | --- | --- | --- |
| Bosonic mode | [`FockSpace`](@extref SecondQuantizedAlgebra.FockSpace) | [`Destroy`](@extref SecondQuantizedAlgebra.Destroy), [`Create`](@extref SecondQuantizedAlgebra.Create) | ``[a, a^\dagger] = 1`` |
| Finite-level system | [`NLevelSpace`](@extref SecondQuantizedAlgebra.NLevelSpace) | [`Transition`](@extref SecondQuantizedAlgebra.Transition) | ``\sigma_{ij}\sigma_{kl} = \delta_{jk}\sigma_{il}`` |
| Qubit in Pauli form | [`PauliSpace`](@extref SecondQuantizedAlgebra.PauliSpace) | [`Pauli`](@extref SecondQuantizedAlgebra.Pauli) | ``\sigma_j\sigma_k = \delta_{jk}I + i\epsilon_{jkl}\sigma_l`` |
| Spin | [`SpinSpace`](@extref SecondQuantizedAlgebra.SpinSpace) | [`Spin`](@extref SecondQuantizedAlgebra.Spin) | ``[S_j,S_k] = i\epsilon_{jkl}S_l`` |
| Canonical quadratures | [`PhaseSpace`](@extref SecondQuantizedAlgebra.PhaseSpace) | [`Position`](@extref SecondQuantizedAlgebra.Position), [`Momentum`](@extref SecondQuantizedAlgebra.Momentum) | ``[p,x] = -i`` |
| Collective finite-level ensemble | [`CollectiveNLevelSpace`](@extref SecondQuantizedAlgebra.CollectiveNLevelSpace) | [`CollectiveTransition`](@extref SecondQuantizedAlgebra.CollectiveTransition) | ``[S^{ij},S^{kl}] = \delta_{jk}S^{il} - \delta_{li}S^{kj}`` |

Use [`ProductSpace`](@extref SecondQuantizedAlgebra.ProductSpace) (written `⊗` or `tensor`) to combine
independent subsystems. For example, a cavity–atom model can contain a [`Destroy`](@extref SecondQuantizedAlgebra.Destroy)
operator acting on subspace `1` and a [`Transition`](@extref SecondQuantizedAlgebra.Transition) operator acting on
subspace `2`. See the [SQA documentation](@extref SecondQuantizedAlgebra :doc:`index`) for the complete operator and
Hilbert-space API.

SQA expressions use ordinary Julia arithmetic. Products are eagerly put into the canonical form
for the selected algebra; `simplify`, [`normal_order`](@extref SecondQuantizedAlgebra.normal_order),
and [`commutator`](@extref SecondQuantizedAlgebra.commutator) provide explicit algebraic operations.

For example, this combines a bosonic cavity mode with a two-level atom using SQA, while Symbolics.jl supplies the scalar parameters:

```julia-repl
julia> using FloquetExpansions

julia> using Symbolics: @variables

julia> cavity = FockSpace(:cavity); atom = NLevelSpace(:atom, 2);

julia> system = cavity ⊗ atom
ℋ(cavity) ⊗ ℋ(atom)

julia> @qnumbers a::Destroy(system, 1) σ::Transition(system, 1, 2, 2);

julia> @variables ωc::Real ωa::Real g::Real A::Real ωd::Real t::Real

julia> H = ωc * a' * a + ωa * σ' * σ + g * (a' * σ + a * σ') + A * cos(ωd * t) * (a + a');

julia> simplify(commutator(a, σ))
0

julia> G = harmonics(H, ωd, t);

julia> G[0]
ωa * σ₂₂ + g * a * σ₂₁ + ωc * a' * a + g * a' * σ₁₂

julia> G[1]
(1//2)*A * a + (1//2)*A * a'
```

Here `H` is an SQA operator expression whose coefficients contain `Symbolics.Num` values.
[`harmonics`](@ref) preserves that exact algebra while converting the time dependence into Floquet
components.

[^complex_symbol]: Use `::Number` rather than `::Complex` when you want one atomic complex parameter;
  `::Complex` represents it as `Complex{Num}` with separate symbolic real and imaginary parts.

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
liouvillian
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

The following constructors encode these channel terms. We allow for both form of a quantum channel: a complete collapse operator or a bare jump operator with a separate rate.

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
rate is nonnegative, nor does constructing a [`liouvillian`](@ref) certify that the result is a
GKLS generator.

### Dissipative frames and Kossakowski coordinates

A [`Liouvillian`](@ref) does not carry a preferred Lindblad representation. When a dissipative
operator frame ``\{F_\mu\}`` is known, FloquetExpansions can extract the coherent Hamiltonian and
Kossakowski matrix exactly in the symbolic operator algebra:

```math
\mathcal L(\rho) = -i[H,\rho]
+ \sum_{\mu\nu} d_{\mu\nu}
\left(F_\mu \rho F_\nu^\dagger
- \frac12\{F_\nu^\dagger F_\mu,\rho\}\right).
```

A [`DissipativeFrame`](@ref) is an ordered set of linearly independent operator directions modulo
the identity. Its ordering fixes the coordinates of the Kossakowski matrix; the frame need not be
orthogonal or normalized. Extraction uses the two-sided sandwich terms directly and verifies the
remaining one-sided action as a Hamiltonian commutator. No Hilbert-space matrix representation or
finite-dimensional cutoff is introduced.

```@docs
DissipativeFrame
kossakowski
hamiltonian
```

The coherent Hamiltonian is defined modulo an additive multiple of the identity, which does not
change its commutator action. Kossakowski coordinates of a raw Floquet expansion use the same frame
and are described in [Floquet expansion](@ref).

### Terms and composition

The left/right representation writes a map as

```math
\mathcal{L}[\rho] = \sum_j c_j A_{j,\mathrm{L}}\rho A_{j,\mathrm{R}}.
```

The [`terms`](@ref) iterator exposes each contribution as a semantic triple
`(left, right, coefficient)`. This gives callers a stable way to inspect a map without depending
on its sparse storage.

```@docs
terms
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
quasienergy calculations and for the Floquet–Liouville analogue in open systems.

```@docs
QuasienergyOperator
```

```@docs
harmonic_range
```

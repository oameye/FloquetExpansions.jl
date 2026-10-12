```@meta
CurrentModule = FloquetExpansions
```

# System

`FloquetExpansions.jl` describes a driven quantum system through its generator. A closed system is
specified by a Hamiltonian whose commutator gives the generator, while an open system is described
directly by a Liouvillian, an SQA `Lindbladian` or `SuperOperator`. The first step for either case is to express the time dependence in a
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

[`floquet_expansion`](@ref) takes the Hamiltonian expression directly. The same action as a map on
density operators is SQA's
`hamiltonian_action`, which returns a
`Lindbladian` without channels.

### Open systems: Lindbladians and superoperators

For an open system, the generator acts directly on density operators:

```math
\frac{d\rho}{dt} = \mathcal{L}(t)[\rho].
```

FloquetExpansions uses the open-system types of SQA rather than a map type of its own. A
`Lindbladian` keeps the structure

```math
\mathcal{L}[\rho]
= -i[H,\rho] + \sum_a \gamma_a\,\mathcal{D}[J_a](\rho),
\qquad
\mathcal{D}[J](\rho)
= J\rho J^\dagger
- \frac{1}{2}\left(J^\dagger J\rho + \rho J^\dagger J\right),
```

and is built from `hamiltonian_action` and
`dissipator` with ordinary arithmetic. A
`SuperOperator` is any polynomial map ``\rho \mapsto \sum_j c_j A_j \rho B_j``; the
`SuperOperator(L)` lowering of a Lindbladian gives its sandwiches. Neither type certifies
positivity: a Lindbladian is a form whose rates may be symbolic or negative, and a general
superoperator need not be of GKSL form. See the
[SQA open-systems manual](https://qojulia.github.io/SecondQuantizedAlgebra.jl/dev/open_systems/) for construction,
composition, and channel normalization.

The two input types differ in what [`floquet_expansion`](@ref) retains. The channels of a
`Lindbladian` are the microscopic input of [`positive_completion`](@ref): they seed the dissipative
frame, and each rate is assumed nonnegative as a whole. Every channel must therefore be diagonal
and unsummed, with a provably real rate, and its Hamiltonian must be Hermitian. A `SuperOperator`
enters the expansion as a bare map, without channels or physical checks.

```julia
using FloquetExpansions
using Symbolics: @variables

h = FockSpace(:cavity); a = Destroy(h, :a)
@variables ω::Real t::Real Δ::Real F::Real κ::Real

L = hamiltonian_action(Δ * a' * a + F * cos(ω * t) * (a + a')) + κ * dissipator(a)
expansion = floquet_expansion(L, ω, t, VanVleck(), 3)
```

The effective components of a Liouvillian expansion are `SuperOperator`s with every finite-level
completeness relation applied, so a component is zero exactly when it is the zero map.

### Dissipative frames and Kossakowski coordinates

A `SuperOperator` does not carry a preferred Lindblad representation. When a dissipative operator
frame ``\{F_\mu\}`` is known, SecondQuantizedAlgebra extracts the coherent Hamiltonian and
Kossakowski matrix exactly in the symbolic operator algebra:

```math
\mathcal L(\rho) = -i[H,\rho]
+ \sum_{\mu\nu} d_{\mu\nu}
\left(F_\mu \rho F_\nu^\dagger
- \frac12\{F_\nu^\dagger F_\mu,\rho\}\right).
```

A `DissipativeFrame` is an ordered set of linearly independent operator directions modulo the
identity. Its ordering fixes the coordinates of the Kossakowski matrix; the frame need not be
orthogonal or normalized. `kossakowski(S, frame)` reads the matrix from the two-sided sandwich
terms and verifies the remaining one-sided action as a Hamiltonian commutator, and
`Lindbladian(S, frame)` returns the map as a `Lindbladian` with one channel per nonzero matrix
entry. No Hilbert-space matrix representation or finite-dimensional cutoff is introduced. These
are SecondQuantizedAlgebra functions, documented in its
[open-systems manual](https://qojulia.github.io/SecondQuantizedAlgebra.jl/dev/open_systems/).

The coherent Hamiltonian `hamiltonian(Lindbladian(S))` is defined modulo an additive multiple of
the identity, which does not change its commutator action. Kossakowski coordinates of a Floquet
expansion use the same frame and are described in [Floquet expansion](@ref).

## Periodic generators

When the generator is periodic with angular frequency ``\omega_d`` and period
``T_d = 2\pi/\omega_d``, write it as

```math
G(t) = \sum_l G_l\,e^{-i l \omega_d t}.
```

Here `l` is an integer harmonic label and ``G_l`` is an operator or `SuperOperator` component. A
[`PeriodicGenerator`](@ref) stores these components together with the drive frequency, and its
component type determines the algebra used by the Floquet expansion. Missing harmonics are zero;
the full time-dependent expression is recovered by evaluating the generator at a symbolic time.

```@docs
PeriodicGenerator
```

### Fourier decomposition

Use [`harmonics`](@ref) to decompose a symbolic Hamiltonian into the convention above. The API
also decomposes a `Lindbladian` or `SuperOperator` by applying Fourier decomposition to the
operators and scalar coefficient of every sandwich, so time dependence may occur in the Hamiltonian, channels, rates, or any
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

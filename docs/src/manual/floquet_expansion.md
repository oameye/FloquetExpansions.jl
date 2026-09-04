```@meta
CurrentModule = FloquetExpansions
CollapsedDocStrings = true
```

# Floquet expansion

Starting from a periodic generator ``\mathcal{G}(t)``, a high-frequency expansion separates
slow evolution from periodic micromotion. If ``\mathcal{V}(t,0)`` denotes the evolution over
time ``t``, the Floquet decomposition is the standard Floquet form
[Shirley1965, Sambe1973](@cite)

```math
\mathcal{V}(t,0) = \mathcal{M}(t)
e^{t\mathcal{G}_\mathrm{eff}}
\mathcal{M}(0)^{-1},
\qquad
\mathcal{M}(t)=e^{\mathcal{K}(t)}.
```

Here ``\mathcal{G}_\mathrm{eff}`` is time independent and ``\mathcal{K}(t)`` is the periodic
generator of the micromotion map. The package currently implements the van Vleck expansion
[VanVleck1929, Eckardt2015](@cite). Hamiltonian and
Liouvillian inputs use the same expansion engine once they have been expressed as
[`PeriodicGenerator`](@ref) values; the corresponding extension to periodic Lindblad generators
is discussed in [Ikeda2021, Schnell2021](@cite).

## Choosing a gauge

The effective generator and micromotion are a gauge-dependent split of the same evolution. The
[`VanVleck`](@ref) gauge fixes the constant part of the micromotion generator by requiring

```math
\langle \mathcal{K}\rangle = 0.
```

This convention makes the effective generator independent of the initial phase used to describe
the drive [Eckardt2015, Bukov2015, Goldman2014](@cite). A stroboscopic Floquet gauge instead chooses a
reference time ``t_0`` and normalizes the two-point micromotion map as
``\mathcal{P}_{t_0}(t)=\mathcal{M}(t)\mathcal{M}(t_0)^{-1}``, so that
``\mathcal{P}_{t_0}(t_0)=I``; its effective generator generally depends on ``t_0``
[Bukov2015](@cite). [`Gauge`](@ref) is the extension point for additional choices.

```@docs
Gauge
```

```@docs
VanVleck
```

## Computing an expansion

Prepare the time dependence as a [`PeriodicGenerator`](@ref), or pass a symbolic Hamiltonian or
Liouvillian together with its drive frequency and time variable. The expansion returns a
[`FloquetExpansion`](@ref) containing the retained effective-generator and micromotion
coefficients together with its finite-realization state.

```@docs
FloquetExpansion
```

```@docs
floquet_expansion
```

## Reading the result

For an expansion truncated at order ``N``, the retained effective generator and micromotion have
the inverse-frequency structure

```math
\mathcal{G}_\mathrm{eff}^{[N]} =
\sum_{n=0}^{N-1}\omega_d^{-n}\mathcal{G}_\mathrm{eff}^{(n)},
\qquad
\mathcal{K}^{[N]} =
\sum_{n=1}^{N-1}\omega_d^{-n}\mathcal{K}^{(n)}.
```

The stored coefficients do not include these powers of the drive frequency. The component
accessors attach the scaling when the result is read. [`effective_component`](@ref) always returns
one retained perturbative contribution. [`effective_generator`](@ref) instead returns the finite
realization: for an [`Uncompleted`](@ref) expansion this is the truncated sum above, while a
positively completed expansion may replace that finite realization without modifying any retained
component.

```@docs
effective_generator
effective_component
```

```@docs
micromotion
```

### Coherent and dissipative components

For a Hamiltonian expansion, [`hamiltonian`](@ref) returns the finite effective Hamiltonian and
[`hamiltonian_component`](@ref) returns one retained inverse-frequency contribution. For a
Liouvillian expansion, the same accessors return its coherent Hamiltonian sector, defined modulo an
additive multiple of the identity.

A raw Liouvillian Floquet expansion does not carry a preferred dissipative representation.
Kossakowski coordinates therefore require an explicit [`DissipativeFrame`](@ref):

```julia
frame = DissipativeFrame(a, a^2)
d = kossakowski(vv, frame)
d1 = kossakowski_component(vv, frame, 1)
```

The component accessors include the corresponding inverse-drive-frequency scaling, consistently
with [`effective_component`](@ref). The operator-frame construction and Liouvillian-level
Kossakowski representation are described in [System](@ref).

```@docs
hamiltonian_component
kossakowski_component
```

## Completion state

Raw Floquet expansions carry [`Uncompleted`](@ref) state. Positive completion is defined as a
new `FloquetExpansion → FloquetExpansion` transformation: the input periodic generator, retained
effective components, micromotion components, gauge, and expansion order remain unchanged, while
the finite no-index effective Liouvillian may be replaced by a completely-positive realization.
This separation keeps the controlled high-frequency data distinct from the higher-order
continuation used to restore complete positivity.

```@docs
Completion
Uncompleted
CompletionAlgorithm
Gram
Spectral
CompletionFactorization
positive_completion
```

`Gram()` and `Spectral()` select the two completion families. This release establishes their type
state and dispatch boundary; the Gram and spectral completion constructions themselves are added in
the corresponding implementation stages. Calling `positive_completion` therefore does not yet
construct a completed generator.

The overload without an explicit dissipative frame is intended for the high-level physical
constructor

```julia
floquet_expansion(H, ωd, t, gauge, order; channels=...)
```

which retains microscopic collapse/jump information internally before lowering to the generic
Liouvillian algebra. An expansion constructed from an already lowered [`Liouvillian`](@ref), from
[`harmonics`](@ref), or directly from a [`PeriodicGenerator`](@ref) does not claim such microscopic
provenance. The explicit-frame overload remains available for those algebraic input paths.

## Order and interpretation

The package uses `order = 1` for the period average. Increasing `order` retains additional
inverse-frequency contributions, but the expansion is asymptotic rather than convergent. Beyond
a problem-dependent optimal order, retaining more terms can make the approximation worse for a
fixed drive.

For dissipative systems, an [`Uncompleted`](@ref) finite-order effective Liouvillian is an
algebraic truncation of the common map expansion and is not guaranteed to retain GKLS form or
complete positivity. Positive completion changes the finite realization without rewriting the
already controlled perturbative components or micromotion.
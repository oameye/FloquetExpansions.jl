```@meta
CurrentModule = FloquetExpansions
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
Liouvillian inputs use the same expansion algorithms once they have been expressed as
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

A truncated van Vleck Liouvillian need not be completely positive. The [`GKSLNormalForm`](@ref)
gauge keeps the oscillatory part of the van Vleck micromotion and adds a static similarity
``W = e^{S}`` with ``S = \sum_n \omega_d^{-n} S_n``. Each ``S_n`` is solved inside the expansion
recurrence so that the dark part of the order-``n`` Kossakowski coefficient is positive
semidefinite. The effective generator ``W^{-1}\mathcal{G}_\mathrm{VV}W`` therefore has the van
Vleck spectrum through the retained order, while every truncation is of GKSL form.

Where every ``S_n`` vanishes the two gauges coincide, which is always the case for Hamiltonian
generators. When the order-``n`` coefficient has a dark part that no static gauge can cancel, a
new channel is born with an exact rate ``\propto \omega_d^{-n}``. The finite generator is assembled
from the graded jump amplitudes, so it agrees with the retained series through the retained order
and differs from it only by a positive remainder beyond it. Bosonic modes are treated without a
Fock cutoff, and all coefficients are exact rationals. With symbolic parameters the expansion
is returned symbolically, with its positivity conditions, when no static gauge is needed through
the requested order, which is then the [`Gram`](@ref) completion of the van Vleck expansion.
An order that needs a static gauge requires numeric parameter values.

```@docs
GKSLNormalForm
```

## Computing an expansion

Prepare the time dependence as a [`PeriodicGenerator`](@ref), or pass a symbolic Hamiltonian or
Liouvillian together with its drive frequency and time variable. The expansion returns a
[`FloquetExpansion`](@ref) containing the effective-generator and micromotion coefficients.

```@docs
FloquetExpansion
order
```

```@docs
floquet_expansion
```

## Reading the result

For an expansion with `order = N`, the retained effective generator and micromotion have the
inverse-frequency structure

```math
\mathcal{G}_\mathrm{eff}^{[N]} =
\sum_{n=0}^{N-1}\omega_d^{-n}\mathcal{G}_\mathrm{eff}^{(n)},
\qquad
\mathcal{K}^{[N]} =
\sum_{n=1}^{N-1}\omega_d^{-n}\mathcal{K}^{(n)}.
```

The distinction between the finite effective generator and its retained perturbative components
becomes important after positive completion; see [Positive completion](@ref positive-completion-manual).

```@docs
effective_generator
effective_component
micromotion
```

### Coherent and dissipative components

For a Hamiltonian expansion, [`hamiltonian`](@ref) returns the finite effective Hamiltonian and
[`hamiltonian_component`](@ref) returns one retained inverse-frequency contribution. For a
Liouvillian expansion, the same accessors return its coherent Hamiltonian sector, defined modulo an
additive multiple of the identity.

A raw Liouvillian Floquet expansion does not carry a preferred dissipative representation.
Kossakowski coordinates therefore require an explicit [`DissipativeFrame`](@ref). The
[`kossakowski_component`](@ref) docstring gives the complete raw-expansion workflow, including the
finite Kossakowski form and an individual retained component.

The operator-frame construction and Liouvillian-level Kossakowski representation are described in
[System](@ref).

```@docs
hamiltonian_component
kossakowski_component
```

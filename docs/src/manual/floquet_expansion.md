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
[`FloquetExpansion`](@ref) containing the effective-generator and micromotion coefficients.

```@docs
FloquetExpansion
```

```@docs
floquet_expansion
```

## Reading the result

For an expansion truncated at order ``N``, the effective generator and micromotion have the
inverse-frequency structure

```math
\mathcal{G}_\mathrm{eff}^{[N]} =
\sum_{n=0}^{N-1}\omega_d^{-n}\mathcal{G}_\mathrm{eff}^{(n)},
\qquad
\mathcal{K}^{[N]} =
\sum_{n=1}^{N-1}\omega_d^{-n}\mathcal{K}^{(n)}.
```

The stored coefficients do not include these powers of the drive frequency. The accessors attach
the scaling when the result is read.

```@docs
effective_generator
```

```@docs
micromotion
```
## Order and interpretation

The package uses `order = 1` for the period average. Increasing `order` retains additional
inverse-frequency contributions, but the expansion is asymptotic rather than convergent. Beyond
a problem-dependent optimal order, retaining more terms can make the approximation worse for a
fixed drive.

For dissipative systems, the finite-order effective Liouvillian is an algebraic truncation of the
common map expansion. It is not projected back into GKLS form, so complete positivity and other
physicality properties require a separate analysis.

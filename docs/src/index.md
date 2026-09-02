```@raw html
---
layout: home

hero:
  name: FloquetExpansions.jl
  text: High-frequency expansions for periodically driven quantum systems
  tagline: Derive effective generators and micromotion symbolically in Julia.
  actions:
    - theme: brand
      text: API Reference
      link: API/
    - theme: alt
      text: Literature
      link: literature/
    - theme: alt
      text: View on GitHub
      link: https://github.com/oameye/FloquetExpansions.jl

features:
  - icon: 〰️
    title: High frequency expansion
    details: Compute a time-independent effective generator together with its periodic micromotion.
  - icon: ∑
    title: Symbolic operator algebra
    details: Build drives from SecondQuantizedAlgebra.jl expressions and retain exact symbolic coefficients throughout the expansion.
  - icon: 🌀
    title: Different gauges
    details: Use a zero-average micromotion generator, making the effective generator independent of the drive's initial phase.
  - icon:
      src: https://sciml.ai/assets/favicon.png
      alt: SciML logo
      wrap: true
    title: Native Julia workflow
    details: Compose the package with Julia's scientific-computing ecosystem including QuantumToolbox.jl, QuantumCumulants.jl and DifferentialEquations.jl.
---
```

```@meta
CurrentModule = FloquetExpansions
```

`FloquetExpansions.jl` computes high-frequency expansions for periodically driven quantum systems. Given a periodic generator ``\mathcal{G}(t)``—Hamiltonian or Liouvillian—it derives a static effective generator ``\mathcal{G}_\mathrm{eff}`` and the generator ``\mathcal{K}(t)`` of a periodic micromotion superoperator ``\mathcal{M}(t)=e^{\mathcal{K}(t)}`` such that

```math
\mathcal{V}(t,0) = \mathcal{M}(t) e^{t\mathcal{L}_\mathrm{eff}} \mathcal{M}(0)^{-1}.
```

The package currently implements the van Vleck expansion, a high-frequency approximation that separates slow effective dynamics from fast micromotion [Eckardt2015](@cite).

## Get Started

Install the package in Julia's package manager:

```julia-repl
pkg> add FloquetExpansions
```

Construct a driven qubit and expand it to the desired order in the inverse drive frequency:

```julia
using FloquetExpansions
using Symbolics: @variables

space = PauliSpace(:qubit)
σz = Pauli(space, :σ, 3)
σx = Pauli(space, :σ, 1)

@variables ω::Real t::Real Δ::Real A::Real γ::Real
H = (Δ / 2) * σz + A * cos(ω * t) * σx

expansion = floquet_expansion(H, ω, t, VanVleck(), 2; channels=(jump(σx, γ),))
L_eff = effective_generator(expansion)
K = micromotion(expansion)(t)
```

See the [API reference](API.md) for the full interface.

## Main Features

### High-frequency expansion

[`floquet_expansion`](@ref) accepts either a periodic generator or a symbolic time-dependent Hamiltonian. Its result exposes [`effective_generator`](@ref) and [`micromotion`](@ref) without materializing a Floquet-Sambe matrix.

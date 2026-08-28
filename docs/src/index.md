```@raw html
---
layout: home

hero:
  name: FloquetExpansions.jl
  text: High-frequency expansions for periodically driven quantum systems
  tagline: Derive effective Hamiltonians and micromotion operators symbolically in Julia.
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
    details: Compute a time-independent effective Hamiltonian together with the kick operator governing micromotion.
  - icon: ∑
    title: Symbolic operator algebra
    details: Build drives from SecondQuantizedAlgebra.jl expressions and retain exact symbolic coefficients throughout the expansion.
  - icon: 🌀
    title: Different gauges
    details: Use a zero-average kick operator, making the effective Hamiltonian independent of the drive's initial phase.
  - icon:
      src: https://sciml.ai/assets/favicon.png
      alt: SciML logo
      wrap: true
    title: Native Julia workflow
    details: Compose the package with Julia's scientific-computing ecosystem including QuantumToolbox.jl, QuantumCumulants.jl and DifferentialEquations.jl.
---
```

- feat: Floquet-Magnus gauge with initial-time condition
- docs: QuantumToolbox and QuantumCumulants interoperability
- design: Floquet-Lindblad effective Liouvillian API and guarantees
- docs: convergence, asymptotic truncation, and high-frequency validity
- research: define scope and reference for flow-equation expansion

```@meta
CurrentModule = FloquetExpansions
```

`FloquetExpansions.jl` computes high-frequency expansions for periodically driven quantum systems. Given a periodic Hamiltonian ``H(t)``, it derives a static effective Hamiltonian ``H_\mathrm{eff}`` and a periodic kick operator ``K(t)`` such that

```math
U(t) = e^{-i K(t)} e^{-i H_\mathrm{eff} t} e^{i K(0)}.
```

The package currently implements the van Vleck expansion, a high-frequency approximation that separates slow effective dynamics from fast micromotion [Eckardt2015](@cite).

## Get Started

Install the package in Julia's package manager:

```julia-repl
pkg> add FloquetExpansions
```

Construct a time-dependent operator and expand it to the desired order in the inverse drive frequency:

```julia
using FloquetExpansions
using SecondQuantizedAlgebra
using Symbolics

h = FockSpace(:cavity)
a = Destroy(h, :a)

@variables w::Real t::Real g::Real
H = \omega * (a' * a) + g * cos(w * t) * (a + a')

expansion = floquet_expansion(H, w, t, VanVleck(), 2)
Heff = effective_hamiltonian(expansion)
K = kick_operator(expansion, t)
```

See the [API reference](API.md) for the full interface.

## Main Features

### High-frequency expansion

[`floquet_expansion`](@ref) accepts either a harmonic representation or a symbolic time-dependent operator. Its result exposes [`effective_hamiltonian`](@ref) and [`kick_operator`](@ref) without materializing a Floquet-Sambe matrix.

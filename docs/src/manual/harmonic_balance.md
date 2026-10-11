```@meta
CurrentModule = FloquetExpansions
```

# [Quantum harmonic balance](@id harmonic-balance-manual)

Harmonic balance describes one driven resonator through several carrier frequencies at once.
First choose the carriers, then build the time-independent carrier generator, then read physical
observables back.

## Choosing the carriers

```@docs
CarrierEmbedding
carrier_modes
```

## Building the carrier generator

```@docs
harmonic_balance
```

## Reading physical observables

```@docs
reconstruct
```

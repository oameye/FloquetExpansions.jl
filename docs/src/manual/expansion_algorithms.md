```@meta
CurrentModule = FloquetExpansions
```

# [Expansion algorithms](@id expansion-algorithms-manual)

After the gauge, the second choice in a Floquet expansion is the algorithm that computes it.

```@docs
ExpansionAlgorithm
```

## Choosing an algorithm

Keep the default unless a high-order expansion is slow; then time [`BlochFeshbach`](@ref) on the
same generator and order.

```@docs
HoriDeprit
BlochFeshbach
```

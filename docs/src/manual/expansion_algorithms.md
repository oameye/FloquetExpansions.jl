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

## Keeping the raw Liouvillian series

For Liouvillian input, both algorithms complete the expansion with [`Gram`](@ref) by default, so
the finite [`effective_generator`](@ref) is a GKLS generator while the retained coefficients and
micromotion stay canonical. Pass `complete_positive=Val(false)` to either algorithm to keep the raw
series, for instance to complete it in a fixed frame or with [`Spectral`](@ref) through
[`positive_completion`](@ref).

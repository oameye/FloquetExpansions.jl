```@meta
CurrentModule = FloquetExpansions
```

# Period map

[`floquet_expansion`](@ref) describes the driven dynamics by an effective generator, which
exists only where the one-period map has a time-independent logarithm. [`period_map`](@ref)
expands the one-period map itself and returns it as a superoperator that is completely positive at
every drive frequency. Use it when complete positivity of the stroboscopic evolution matters more than a
generator, or when no GKSL effective generator exists.

```@docs
period_map
PeriodMap
```

## Reading the result

The superoperator is the map. The trace defect measures how far it is from trace
preserving; it starts beyond the retained order.

```@docs
superoperator
trace_defect
order(::PeriodMap)
```

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

## Kraus form

The same map is a sum over jump records of squared amplitudes, integrated over the jump times.
Each amplitude is an operator whose coefficients depend on the fractional jump times, so the map
is in continuous Kraus form and complete positivity is manifest.

```@docs
jump_amplitudes
JumpAmplitude
weight
amplitude
jump_times
```

## Representation changes

Mixing jump channels of equal rate by a unitary leaves the map unchanged. An affine shift
``L \to L + c`` with the compensating Hamiltonian describes the same dynamics but changes which
histories count as jumps, so the two period maps agree only through the retained order.

```@meta
CurrentModule = FloquetExpansions
CollapsedDocStrings = true
```

# Positive completion

A finite-order high-frequency expansion of a periodic Liouvillian is an algebraic truncation of
the exact Floquet generator. Even when the microscopic evolution is completely positive, the
finite truncated effective Liouvillian need not itself retain GKLS form. Positive completion is
the separate finite-realization layer that restores a completely positive effective generator
without rewriting the perturbative data already fixed by the expansion.

For an expansion truncated at order ``N``, the completion is constructed so that the completed
finite generator agrees with the raw high-frequency expansion through the retained order. In
particular, the retained effective components and micromotion remain the controlled perturbative
data; completion supplies only a higher-order continuation of the finite effective Liouvillian.

## Completion state

Raw Floquet expansions carry [`Uncompleted`](@ref) state. Positive completion is represented as a
new `FloquetExpansion → FloquetExpansion` transformation. The input periodic generator, retained
effective components, micromotion components, gauge, and expansion order are preserved, while the
finite no-index effective Liouvillian may be replaced by a completely positive realization.

This distinction is reflected directly in the public accessors:

- [`effective_component`](@ref) always returns the original retained perturbative coefficient;
- [`micromotion`](@ref) always returns the original retained micromotion;
- [`effective_generator`](@ref) returns the current finite realization, raw or completed.

Thus positive completion changes the finite model, not the already-controlled Floquet data.

```@docs
Completion
Uncompleted
CompletionAlgorithm
CompletionFactorization
positive_completion
```

## Completion algorithms

[`Gram`](@ref) and [`Spectral`](@ref) are the public algorithm selectors for the two completion
families.

```@docs
Gram
Spectral
```

This implementation stage establishes the completion state, public algorithm types, and dispatch
boundary. The Gram and spectral completion constructions themselves are implemented in the
corresponding follow-up stages. Until then, `positive_completion` establishes the correct public
entry points but does not yet construct a completed generator.

## Physical and algebraic input paths

Automatic completion requires enough microscopic information to identify the dissipative sector.
The high-level physical constructor

```julia
floquet_expansion(H, ωd, t, gauge, order; channels=...)
```

therefore retains the ordered collapse/jump input internally before lowering the dynamics to the
generic Liouvillian algebra. This provenance belongs to the high-level construction layer only;
[`Liouvillian`](@ref) and [`PeriodicGenerator`](@ref) remain basis-free algebraic objects.

An expansion constructed from an already lowered [`Liouvillian`](@ref), from
[`harmonics`](@ref), or directly from a [`PeriodicGenerator`](@ref) does not claim microscopic
channel provenance. For those paths, completion requires an explicit [`DissipativeFrame`](@ref).

The frame itself is an ordered operator basis for the dissipative module, as described in
[System](@ref). The physical collapse/jump semantics, including real time-periodic rates, are also
documented there rather than in the completion layer.

## Micromotion

Positive completion does not require a modification of the retained kick operator. If the raw
expansion is controlled through order ``N`` and completion changes the finite effective generator
only beyond that order, then the same retained micromotion remains valid through the claimed
order. The completed approximation therefore uses the original micromotion together with the
completed finite effective generator.

A change of micromotion would only be required if one demanded exact equality with the original
finite truncated propagator at all times rather than perturbative agreement through the retained
order.

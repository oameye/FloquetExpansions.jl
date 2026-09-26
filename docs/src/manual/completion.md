```@meta
CurrentModule = FloquetExpansions
```

# [Positive completion](@id positive-completion-manual)

A finite-order high-frequency expansion of a periodic Lindblad generator need not itself be a
GKLS generator, even when the microscopic dynamics is Markovian and completely positive
[Schnell2021, Ikeda2021, Haddadfarshi2015](@cite). Positive completion supplies a finite completely positive continuation
without changing the retained Floquet coefficients or micromotion. It only fixes dissipative information beyond
the retained order. Different completion methods can therefore agree with the same retained
expansion while differing at higher order.

For Liouvillian Van Vleck expansions, [`HoriDeprit`](@ref) and [`BlochFeshbach`](@ref) use
`Gram()` completion by default. Request the raw canonical HFE with
`complete_positive=Val(false)`. The explicit [`positive_completion`](@ref) API remains useful for
fixed dissipative frames, `Spectral()`, and completion-method comparisons.

```@docs
positive_completion
```

## Choosing a completion algorithm

The package provides two CP completion algorithms: [`Gram`](@ref) and [`Spectral`](@ref).

Both methods preserve the retained Kossakowski coefficients, but their finite positive
continuations need not agree beyond the retained order.

```@docs
Gram
Spectral
```

## Reading a completed expansion

A completed result remains a [`FloquetExpansion`](@ref) struct. The [`effective_component`](@ref),
[`hamiltonian_component`](@ref), and [`micromotion`](@ref) expose retained Floquet data. Whereas,
[`effective_generator`](@ref), [`hamiltonian`](@ref), [`kossakowski`](@ref), and [`channels`](@ref)
expose the finite completed model.

```@docs
channels
dissipative_frame
```

For a raw expansion, Kossakowski coordinates require an explicit [`DissipativeFrame`](@ref). A
completed expansion stores the finalized frame, so the no-frame forms of [`kossakowski`](@ref) and
[`kossakowski_component`](@ref) are available directly.

## Symbolic parameter conditions

Symbolic completion can depend on parameters whose sign or nonzero status is not decidable
structurally. The two kinds of assumptions are kept separate:

- [`positivity_conditions`](@ref) records physical nonnegativity assumptions;
- [`regularity_conditions`](@ref) records nonzero assumptions defining the current fixed-rank
  symbolic stratum.

A regularity assumption is not a positivity assumption. When a required pivot vanishes, the
factorization belongs to a different rank stratum and must be reconsidered there.

```@docs
positivity_conditions
regularity_conditions
```

## Factorization diagnostics

The common physical API is algorithm independent. Method-specific diagnostics are
available through [`factorization`](@ref): Gram completion records graded amplitudes and onset
history, while spectral completion records completed rates, perturbative branch vectors, and
Puiseux onset metadata.

```@docs
factorization
GramFactorization
SpectralFactorization
```

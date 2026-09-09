```@meta
CurrentModule = FloquetExpansions
CollapsedDocStrings = true
```

# [Positive completion](@id positive-completion-manual)

A finite-order high-frequency expansion of a periodic Lindblad generator need not itself be a
GKLS generator, even when the microscopic dynamics is Markovian and completely positive
[Schnell2021](@cite). Positive completion supplies a finite completely positive continuation
without changing the retained Floquet coefficients or micromotion.

The normal workflow is explicit:

```julia
vv = floquet_expansion(H, ω, t, VanVleck(), order; channels=(jump(J, γ),))
cp = positive_completion(vv, Gram())
```

`floquet_expansion` remains raw by default. Completion only fixes dissipative information beyond
the retained order. Different completion methods can therefore agree with the same retained
expansion while differing at higher order.

```@docs
positive_completion
```

The mathematical construction and its limitations are developed in
[CP-preserving completion](@ref cp-preserving-completion-theory). The driven-qubit and bosonic
calculations are worked out in [CP-completion examples](@ref cp-completion-examples).

## Choosing a completion algorithm

Use [`Gram`](@ref) for the general algebraic construction and [`Spectral`](@ref) when a suitable
leading spectral frame is available. Supplying a [`DissipativeFrame`](@ref) explicitly is useful
when comparing algorithms in the same physical operator basis.

```julia
cp_gram = positive_completion(vv, Gram(), frame)
cp_spectral = positive_completion(vv, Spectral(), spectral_frame)
```

Both methods preserve the retained Kossakowski coefficients, but their finite positive
continuations need not agree beyond the retained order.

```@docs
Gram
Spectral
```

## Reading a completed expansion

A completed result remains a [`FloquetExpansion`](@ref). The retained perturbative accessors retain
the same meaning as before completion:

```julia
effective_component(cp, n)
hamiltonian_component(cp, n)
micromotion(cp)
```

The finite accessors describe the completed model:

```julia
effective_generator(cp)
hamiltonian(cp)
kossakowski(cp)
channels(cp)
```

The completed channels reconstruct the finite generator together with its coherent sector,

```julia
liouvillian(hamiltonian(cp); channels=channels(cp)) == effective_generator(cp)
```

```@docs
channels
dissipative_frame
```

For a raw expansion, Kossakowski coordinates require an explicit [`DissipativeFrame`](@ref). A
completed expansion stores the finalized frame, so the no-frame forms of `kossakowski` and
`kossakowski_component` are available directly:

```julia
frame = DissipativeFrame(a, a^2)
d_raw = kossakowski(vv, frame)

cp = positive_completion(vv, Gram(), frame)
d_cp = kossakowski(cp)
d1 = kossakowski_component(cp, 1)
```

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

The common physical API is algorithm independent. Method-specific scientific diagnostics are
available through [`factorization`](@ref): Gram completion records graded amplitudes and onset
history, while spectral completion records completed rates, perturbative branch vectors, and
Puiseux onset metadata.

```@docs
factorization
GramFactorization
SpectralFactorization
CompletionObstruction
FractionalJumpOnset
```

These diagnostic types are stable qualified API rather than exported names. For example, use
`FloquetExpansions.GramFactorization` when dispatching on the returned factorization type.

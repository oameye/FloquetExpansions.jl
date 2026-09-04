```@meta
CurrentModule = FloquetExpansions
CollapsedDocStrings = true
```

# Positive completion

A finite-order high-frequency expansion of a periodic Lindblad generator need not itself produce
a GKLS effective generator, even when the microscopic dynamics is Markovian and completely
positive [Schnell2021](@cite). Positive completion restores complete positivity without changing
the perturbative information fixed by the high-frequency expansion.

For `order = N`, write the retained Kossakowski series as

```math
d^{[N]}(\epsilon)=\sum_{n=0}^{N-1}\epsilon^n d^{(n)},
\qquad \epsilon=\omega_d^{-1}.
```

A positive completion constructs a finite matrix ``\widetilde d^{[N]}(\epsilon)`` such that

```math
\widetilde d^{[N]}(\epsilon)\succeq0,
\qquad
\widetilde d^{[N]}(\epsilon)-d^{[N]}(\epsilon)=\mathcal O(\epsilon^N).
```

Thus the completion modifies only terms beyond the retained order. The perturbative effective
components remain unchanged, while [`effective_generator`](@ref) returns the completed finite
generator.

## Micromotion

The retained kick operator does not need to be modified. If
``\mathcal L_{\rm CP}^{[N]}-\mathcal L_{\rm eff}^{[N]}=\mathcal O(\epsilon^N)``, then the same
retained micromotion gives a completed approximation with the same accuracy,

```math
\Phi_{\rm CP}^{[N]}(t)
=
e^{\mathcal K^{[N]}(t)}
e^{t\mathcal L_{\rm CP}^{[N]}}
e^{-\mathcal K^{[N]}(0)}.
```

A change of micromotion is only required if one asks for exact equality with the original finite
truncated propagator rather than agreement through the retained order.

## API

```@docs
Completion
Uncompleted
CompletionAlgorithm
Gram
Spectral
CompletionFactorization
positive_completion
```

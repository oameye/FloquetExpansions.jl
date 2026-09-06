```@meta
CurrentModule = FloquetExpansions
CollapsedDocStrings = true
```

# Positive completion

A finite-order high-frequency expansion of a periodic Lindblad generator need not itself be a
GKLS generator, even when the microscopic dynamics is Markovian and completely positive
[Schnell2021](@cite). Positive completion supplies a finite completely-positive continuation
without changing the perturbative information that the Floquet expansion has already fixed.

The normal workflow is explicit:

```julia
vv = floquet_expansion(H, ω, t, VanVleck(), order; channels=(jump(J, γ),))
cp = positive_completion(vv, Gram())
```

`floquet_expansion` remains raw by default. Completion does not rewrite the retained effective
components, coherent components, or micromotion. It only chooses higher-order dissipative data that
were not fixed by the truncation.

This distinction is important. A finite positive completion is generally **not unique** beyond the
retained order, and complete positivity of the completed finite model is not a convergence theorem
for the high-frequency series or a statement that the exact Floquet logarithm is GKLS. The
completed model is a controlled higher-order continuation of the retained perturbative data.

For a truncation through order ``N`` (`order = N + 1` in the API), write the retained Kossakowski
series in a chosen dissipative frame as

```math
d^{[N]}(\epsilon)=\sum_{n=0}^{N}\epsilon^n d^{(n)},
\qquad \epsilon=\omega_d^{-1}.
```

Here ``\omega_d`` is the positive physical angular drive frequency. This orientation convention is
irrelevant for even inverse-frequency powers but matters when a positive dissipative rate first
appears at odd order. In that case `Spectral()` records ``\omega_d`` in both the positivity and
regularity condition sets, so the strict assumption ``\omega_d>0`` is explicit rather than hidden.

A completion constructs a finite Hermitian form ``\widetilde d^{[N]}`` such that

```math
\widetilde d^{[N]}\succeq0,
\qquad
\widetilde d^{[N]}-d^{[N]}=\mathcal O(\epsilon^{N+1}).
```

See [CP-completion examples](@ref) for the driven-qubit and bosonic calculations used to validate
these statements analytically.

## Reading a completed expansion

A completed result uses the same [`FloquetExpansion`](@ref) interface as a raw result:

```julia
effective_generator(cp)
effective_component(cp, n)
hamiltonian(cp)
hamiltonian_component(cp, n)
kossakowski(cp)
kossakowski_component(cp, n)
dissipative_frame(cp)
channels(cp)
positivity_conditions(cp)
regularity_conditions(cp)
factorization(cp)
```

The finite generator returned by [`effective_generator`](@ref) is the completed generator, whereas
[`effective_component`](@ref) continues to return the original retained Van Vleck coefficients.
Likewise, [`hamiltonian_component`](@ref) and [`micromotion`](@ref) are retained perturbative data.

For a raw expansion, Kossakowski coordinates require an explicit [`DissipativeFrame`](@ref):

```julia
frame = DissipativeFrame(a, a^2)
d_raw = kossakowski(vv, frame)
d1 = kossakowski_component(vv, frame, 1)
```

A completed expansion stores its finalized frame, so the explicit frame is no longer required:

```julia
cp = positive_completion(vv, Gram(), frame)
d_cp = kossakowski(cp)
```

The completed result owns this finalized representation rather than retaining mutable caller-owned
frame data. It also caches the retained physical Kossakowski coefficients and coherent Hamiltonian
used to build the completed generator, so completed accessors do not repeat GKSL coordinate
extraction. Public representation accessors remain defensive against mutation.

The two-argument form `positive_completion(vv, algorithm)` performs symbolic automatic frame
discovery and is the convenient frontend. Because the number of independent generated directions
is determined algebraically at runtime, its frame arity is dynamic. When a frame is known, the
three-argument form `positive_completion(vv, algorithm, frame)` is the inference-oriented
computational core and should be preferred in performance-sensitive code.

The common physical reconstruction contract is

```julia
liouvillian(hamiltonian(cp); channels=channels(cp)) == effective_generator(cp)
```

up to the package's exact symbolic simplification conventions. The returned channels therefore
represent the finite completed generator, not merely the retained truncation.

## Choosing a completion algorithm

[`Gram`](@ref) is the general algebraic path. [`Spectral`](@ref) is a restricted decay-rate/HCM
realization that is useful when a suitable spectral frame is already known.

```julia
cp_gram = positive_completion(vv, Gram())
cp_gram_frame = positive_completion(vv, Gram(), frame)

cp_spectral = positive_completion(vv, Spectral())
cp_spectral_frame = positive_completion(vv, Spectral(), frame)
```

The two methods must agree with the raw expansion through the retained order, but their finite
higher-order positive continuations may differ.

### Gram completion

[`Gram`](@ref) constructs a graded collapse-amplitude matrix

```math
B^{[N]}(\epsilon)=\sum_p\epsilon^p B^{(p)}
```

such that

```math
\Pi_N\!\left(B^{[N]}B^{[N]\dagger}\right)=d^{[N]}.
```

The completed Kossakowski matrix is the **untruncated** product

```math
\widetilde d^{[N]}=B^{[N]}B^{[N]\dagger}\succeq0.
```

The leading Hermitian form is split algebraically into active and dark sectors. The active sector is
factorized by graded ``LDL^\dagger`` elimination and scalar series square roots. Coupling to dark
directions is handled through solves and a Feshbach/Schur residual. No symbolic Kossakowski
eigendecomposition, Hilbert-space matrix representation, or Liouville-space vectorization is
required.

This is why a non-diagonal Cartesian Pauli frame is a perfectly valid input to `Gram()`. The method
does not need a spectral frame first.

When microscopic [`collapse`](@ref) or [`jump`](@ref) channels are available, automatic frame
construction preserves their order and appends independent Floquet-generated dissipative
directions as needed. An explicit [`DissipativeFrame`](@ref) fixes the representation when the
scientific basis itself matters.

[`factorization`](@ref) returns a [`GramFactorization`](@ref). Its graded `amplitudes` encode a jump
amplitude gauge for the completed Gram form, while `onsets` records when independent dissipative
channels first open.

### Spectral/HCM completion

[`Spectral`](@ref) works in a restricted frame in which the leading Kossakowski form is diagonal.
Perturbative Rayleigh--Schrödinger recursion follows decay-rate branches between distinct leading
sectors, and each retained scalar rate is completed by the HCM square construction. Mixing inside
an unresolved degenerate leading sector is rejected rather than hidden behind symbolic
diagonalization.

The completed form has the branch representation

```math
\widetilde d=\sum_a\widetilde\lambda_a\,\phi_a\phi_a^\dagger,
\qquad \widetilde\lambda_a\ge0.
```

A positive branch may begin at an odd power of ``\epsilon``. Its physical rate remains an ordinary
integer-power expression, so [`channels`](@ref) remains valid, although folding the rate into a
collapse amplitude would require a Puiseux power. [`SpectralFactorization`](@ref) records this with
its onset and `puiseux` metadata. For an odd onset, positivity of the finite HCM rate also uses the
physical orientation ``\omega_d>0``; the symbolic implementation records ``\omega_d`` in both
[`positivity_conditions`](@ref) and [`regularity_conditions`](@ref).

Use `factorization(cp)` when algorithm-specific diagnostics are needed. The common physical API is
otherwise deliberately algorithm independent.

## Positivity and regularity conditions

Symbolic completion can depend on parameters whose sign or nonzero status is not decidable
structurally. These assumptions are kept separate:

- [`positivity_conditions`](@ref) contains physical inequalities such as ``\gamma\ge0`` required
  for a nonnegative Kossakowski form;
- [`regularity_conditions`](@ref) contains nonzero assumptions needed to remain on the selected
  fixed-rank symbolic stratum.

A regularity condition is not a positivity assumption. At a parameter value where a required pivot
vanishes, the dissipative rank changes and the completion must be reconsidered on the new stratum.
When the same symbol occurs in both sets, both statements are required. In particular, an odd
spectral rate onset records ``\omega_d\ge0`` and ``\omega_d\ne0`` separately; together these are the
strict physical convention ``\omega_d>0``.

`jump(J, γ)` uses ``\gamma`` as a physical rate: it must be real and nonnegative on the declared
parameter domain.

## Micromotion is retained

Suppose

```math
\mathcal L_{\rm CP}^{[N]}-\mathcal L_{\rm eff}^{[N]}
=\mathcal O(\epsilon^{N+1}).
```

Then no retained kick correction is required. The same micromotion gives

```math
\Phi_{\rm CP}^{[N]}(t)
=
e^{\mathcal K^{[N]}(t)}
e^{t\mathcal L_{\rm CP}^{[N]}}
e^{-\mathcal K^{[N]}(0)},
```

with the same accuracy through order ``N``. A modified kick would only be required if one demanded
exact equality with the original finite truncated propagator rather than perturbative agreement
through the claimed order.

The active/dark construction, onset filtration, frame covariance, and the relation between Gram and
spectral gauges are developed in
[High-frequency expansion](../theory/high_frequency_expansion.md). The exact Floquet interpretation
and its limitations are discussed in [Floquet theory](@ref).

## Development references

The implementation roadmap is tracked in
[FloquetExpansions.jl#54](https://github.com/oameye/FloquetExpansions.jl/issues/54), with the
architecture and terminology recorded in
[PR #55](https://github.com/oameye/FloquetExpansions.jl/pull/55). The analytical and cross-method
validation is in [PR #77](https://github.com/oameye/FloquetExpansions.jl/pull/77); that validation
also exposed the Liouvillian Lie-transform convention corrected under
[issue #79](https://github.com/oameye/FloquetExpansions.jl/issues/79).

## API

```@docs
Completion
Uncompleted
CompletionAlgorithm
Gram
Spectral
CompletionFactorization
CompletionObstruction
FractionalJumpOnset
GramFactorization
SpectralFactorization
positive_completion
dissipative_frame
channels
positivity_conditions
regularity_conditions
factorization
```

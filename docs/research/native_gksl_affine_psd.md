# Native-GKSL affine PSD static section

This research note records the finite-dimensional retained-block problem used by #368.

At generator order `n`, after the ordinary Floquet harmonic elimination and subtraction of already-known Gram products, the dark quotient must satisfy

```math
P_n = \Delta_n + \Phi_n(S_n) \succeq 0.
```

The previously certified implementation chose the minimum-norm coordinates

```math
S_n^{(0)} = -\Phi_n^\#\Delta_n
```

and accepted the order only when the corresponding residual was already positive semidefinite. That is sufficient but not necessary. The actual local solvability criterion is

```math
(\Delta_n + \operatorname{Ran}\Phi_n) \cap \mathrm{PSD}(D_n) \neq \varnothing.
```

The internal `positive_affine_section` helper preserves the old minimum-norm branch exactly when it is already PSD. Otherwise it computes the minimum-Frobenius-norm point in the affine/PSD intersection by Dykstra projection between the Hermitian affine slice and the PSD cone. Hermitian matrix coordinates are chosen as a Frobenius-isometric real vectorization, so the affine projection is an ordinary Euclidean least-squares projection.

Once a feasible residual `P_n` is found, newborn jump amplitudes are still obtained by the shared `positive_gram_factor` kernel. No post-hoc positive-completion workflow is involved.

This is research/internal infrastructure only. It does not yet constitute the full production native-GKSL static solver; BF/HD recurrence integration and physical qutrit/Kerr certification remain follow-up steps.

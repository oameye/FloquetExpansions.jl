# Native-GKSL resolvent section on charged sectors

This research note records the static-slot section used by `test/gksl/native_kerr_graded_slot.jl`. It connects the localized Kerr inverse of #376/#377 to the generic native BF/HD recurrence.

## Problem

The generic slot solves the dark equation

```math
\Delta_n + \Phi_n(S_n) = P_n \succeq 0,
\qquad
\Phi_n(S) = (\pi_n \otimes \bar\pi_n)\, c([\bar{\mathcal L}, S]),
```

and picks the minimum-norm point of the affine PSD slice. The minimum-norm rule depends on a Frobenius metric on the gauge algebra. That metric does not exist in a localized or symbolic coefficient algebra, so an algebraic inverse such as the Kerr ladder `(D+J)^{-1}` cannot reproduce it. Plugging the ladder into the minimum-norm slot is therefore not well defined.

## Grading

Let `Q = (out_L - in_L) - (out_R - in_R)` be the superoperator charge. For a U(1)-symmetric averaged generator, `ad_Lbar` preserves `Q`. The Kossakowski map `c` and the dark projector `\pi_n` also preserve it as long as the leading channels span a graded subspace. This holds inductively, because the newborn form `P_n` is neutral (see below). Hence `\Phi_n` is block diagonal in `Q`. The charged blocks `Q != 0` carry no PSD constraint: pinching a PSD form onto its neutral block keeps it PSD and lowers its norm. So every feasible point can be moved to one with `P_n^{(Q != 0)} = 0`.

## Section

On the regular charged sectors we fix `S_n` by the resolvent condition

```math
[\bar{\mathcal L}, S_n^{\rm ch}] = -\Pi_{Q\neq0}\, \mathcal D[\pi_n \Delta_n \pi_n^\dagger].
```

This removes exactly the charged dark-dark target and adds nothing to the bright or Hamiltonian parts of `E_n`. The neutral block keeps the dense affine PSD section over the neutral gauge algebra (16 of the 72 gauge directions at Fock cutoff 3). The section is defined without any metric. It is intrinsic, so BF and HD reproduce it independently. With no jumps `\Delta_n = 0`, so `S_n = 0` and the Hamiltonian branch is untouched.

For driven Kerr with loss every charged Bohr multiplier is nonzero on the physical lattice. The real part is `-\kappa(q_L+q_R)/2`, and on the neutral-real levels `q_L = -q_R = q` the multiplier is `-2iq[\Delta + \chi(n_L+n_R)]`. A physical zero of the latter is a genuine resonance, and the ladder rejects it. Away from resonance, `ad_Lbar` is invertible on the charged sectors and the triangular ladder `s_t = (y_t - J s_{t+1})/d_t` is its exact inverse.

## Certified on finite-Fock driven Kerr

- The localized ladder agrees with a dense LU of `ad_Lbar` restricted to `Q != 0`, both as an inverse and inside the full recurrence through order 4, at cutoffs 3 and 4.
- BF equals the independent HD recurrence through order 4, together with prefix stability, coefficientwise Gram reconstruction, TP and finite-truncation GKSL positivity.
- No charged dark-dark target survives in any `E_n`.
- Through order 3 the charged target vanishes and the section equals the minimum-norm one. At order 4 the two sections share the neutral gauge, the newborn Gram form and the dark-dark coefficient. They differ by `\delta S \in \ker\Phi_4`, with `E_4 - E_4^{\rm mn} = [\bar{\mathcal L}, \delta S]`, and the spectra of their truncated generators agree to `O(\varepsilon^5)` or better.
- The Hamiltonian branch (`\kappa = 0`) reproduces the ordinary result.

The section is therefore a different but orbit-equivalent GKSL representative. It replaces "minimum-norm on the charged blocks" in Theorem 12, item 3 of the derivation note. That note still states the minimum-norm convention and should be updated if this section is adopted.

This is research test infrastructure only. No public API is added and positive completion is never called.

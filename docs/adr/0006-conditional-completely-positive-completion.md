# Explicit, conditional completely-positive completion

Finite-order Floquet Liouvillian expansions are returned algebraically by default and are not implicitly repaired into GKSL form. Complete positivity is an explicit opt-in post-processing step applied to a `FloquetExpansion`.

`positive_completion(vv, algorithm)` constructs a finite positive continuation that preserves every retained Floquet coefficient through the requested order. The completion may depend on symbolic positivity and regularity conditions. These conditions are reported rather than guessed when the package cannot prove them from microscopic rate provenance or structural algebra.

The first supported completion algorithms are an algebraic Gram/Feshbach construction and a restricted perturbative spectral/HCM construction. Completion does not modify the retained micromotion or Hamiltonian coefficients through the controlled order, and it does not claim that the finite completed model is the exact Floquet GKSL logarithm.

This ADR supersedes the earlier decision that positivity completion was only a future possibility; the feature remains explicit rather than automatic.

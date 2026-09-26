# Default, conditional completely-positive completion

Finite-order Floquet Liouvillian expansions may leave the GKSL cone even when the microscopic driven model is Lindbladian. Van Vleck open-system expansions therefore use graded `Gram()` positive completion by default. The selected Hori–Deprit or Bloch/Feshbach backend still constructs the same canonical retained Floquet coefficients and micromotion; completion changes only the finite effective-generator realization beyond the retained order.

The raw canonical expansion remains available explicitly through `HoriDeprit(; complete_positive=Val(false))` and `BlochFeshbach(; complete_positive=Val(false))`. This is the path for algebraic diagnostics, perturbative comparison, custom completion, and validation against the historical raw HFE.

`positive_completion(vv, algorithm)` remains a public operation on an explicitly raw expansion. It is required for selecting `Spectral()`, supplying a fixed `DissipativeFrame`, or comparing completion algorithms. Completion may depend on symbolic positivity and regularity conditions; unresolved conditions are reported rather than guessed.

The default `Gram()` continuation preserves every retained Floquet coefficient and the retained micromotion. It does not claim that the finite completed generator is the exact Floquet GKSL logarithm. Native CK/open-leg reconstruction of a finite CPTP period map is a separate construction and is not identified with this static completion.

## Gate

`make test` runs the testsets that hold this decision:

- `test/expansion_algorithm_cp_policy.jl`: default selectors equal explicit CP completion and explicit `Val(false)` preserves the raw HD/BF HFE;
- `test/gram_completion.jl`: symbolic positivity and regularity conditions remain distinct;
- `test/cp_completion_validation.jl`: spectral symbolic rank strata separate positivity and regularity.

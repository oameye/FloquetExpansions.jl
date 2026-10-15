# Open-system generators and their GKSL coordinates belong to SecondQuantizedAlgebra

SecondQuantizedAlgebra now owns a symbolic open-system layer: `Lindbladian` keeps `-i[H,·] + Σ γ 𝒟[J, K]` with its channels, and `SuperOperator` is any collected sum of sandwiches `ρ ↦ c AρB`. The second is the representation of ADR 0007, and SQA also owns the multiplication, normal ordering and coefficient algebra beneath it. FloquetExpansions therefore drops its own `Liouvillian` map type and the constructors `hamiltonian_action`, `dissipator`, `liouvillian`, `collapse` and `jump`, together with `compose` and `terms`. It uses the SQA types and names instead. Two map algebras with the same names would collide under `@reexport using SecondQuantizedAlgebra` and would drift apart.

## Representation

`PeriodicGenerator{SuperOperator}` is the Liouvillian component type. `harmonics` lowers a `Lindbladian` to its `SuperOperator` and decomposes every sandwich side and coefficient, as before. SQA keeps the ground projector of an `NLevelSpace` atomic, but the canonical-key amendment of ADR 0007 needs map equality. FloquetExpansions therefore applies `expand_completeness` to every map it produces: harmonics, the composition and commutator of the expansion engine, lifted native results, and completed generators. A component is then zero exactly when it is the zero map. Internally the sandwiches are read as `(left, right, coefficient)` through SQA's `left_operator` and `right_operator`, and maps are built with `SuperBuilder`.

## GKSL coordinates

Reading a map back in GKSL form is a property of the map, not of a Floquet expansion, so it lives in SQA too. SQA owns `JumpBasis`, `kossakowski(S, basis)`, `Lindbladian(S[, basis])`, `support_basis`, `independent_basis` and `GKSLFormError`, ported from FloquetExpansions without change of algorithm. FloquetExpansions' former `DissipativeFrame` is SQA's `JumpBasis`, and its accessor is renamed `jump_basis` so that one name serves both packages. The retained coherent sector is read as `hamiltonian(Lindbladian(S, basis))`. FloquetExpansions keeps only the coefficient-matrix helpers its completion layer needs, and builds every completed model with SQA's `liouvillian(H, J; rates)`. An empty basis is allowed, since it is the basis of a coherent map and automatic discovery already produced one. A `JumpBasis` iterates and indexes like a vector of its operators and has no mutating API, so completion stores the caller's basis without copying it.

SQA rejects floating-point maps in exact extraction unless an `rtol` is passed. FloquetExpansions passes none, so a float model is rejected at completion rather than read up to an unstated tolerance.

## Public introspection

The native lowering and the channel-rate check read SQA values only through its public API: `coefficient_terms` and `has_coefficient_terms` give the exact or floating scalar and the atom powers of a coefficient, and `n_levels`, `ground_state`, `transition_levels`, `operator_axis`, `operator_name` and `acts_on` describe an operator. Lifted operators are rebuilt with the documented constructors and `set_acts_on`. FloquetExpansions previously used SQA's coefficient tails, scalar packing and positional operator constructors, which an internal refactor of SQA could break without a breaking release.

## Microscopic provenance

The physical input of positive completion is a `Lindbladian`. It replaces the `H; channels=(collapse(…), jump(…))` form of ADR 0003. `floquet_expansion(L::Lindbladian, ωd, t, gauge, order)` and `period_map(L::Lindbladian, ωd, t)` read `channels(L)` in SQA's deterministic `channel_order_key` order, which now orders the basis seeds instead of user order. Each rate is checked as `jump` checked it: it must be provably real, a negative numeric rate is rejected, and a symbolic rate is asserted nonnegative as a whole. A correlated channel `𝒟[J, K]` or a channel summed over an index has no single microscopic rate and is rejected with a hint to pass `SuperOperator(L)`. The Hamiltonian of `L` must be Hermitian. The collapse/jump distinction disappears, because SQA normalizes a channel by moving the numeric content of its jump into the rate. A `SuperOperator` input enters without provenance, as an explicitly constructed `Liouvillian` did.

## Completed models

A positively completed or `GKSLNormalForm` expansion stores its finite model as a `Lindbladian`, returned by `Lindbladian(expansion)`, whose `SuperOperator` lowering is `effective_generator(expansion)`. `channels(expansion)` returns its `channel => rate` pairs in `channel_order_key` order. The order is deterministic, which SQA's `Dict` is not. `hamiltonian`, `kossakowski` and `channels` are the SQA functions. FloquetExpansions adds methods only on `FloquetExpansion`, which it owns. The Hamiltonian of a bare map is `hamiltonian(Lindbladian(S))`; a `hamiltonian(S::SuperOperator)` method in FloquetExpansions would be type piracy.

## Considered options

- **Keep `collapse`/`jump` as the provenance API and build an SQA `Lindbladian` internally.** Rejected by the maintainer: it keeps a second channel vocabulary beside SQA's, and FE's `jump` would shadow SQA's public `jump` accessor.
- **Return SQA's channel `Dict` from `channels(expansion)`.** Rejected: its iteration order depends on hashing, so printed results and positional access are not reproducible. SQA's own `channels(L)` now returns the sorted vector.
- **Keep GKSL coordinates in FloquetExpansions.** Rejected: SQA exports `Lindbladian` and `SuperOperator`, and the inverse map from a `SuperOperator` to a `Lindbladian` belongs beside them. Two `kossakowski` functions would collide under `@reexport`.

## Gate

`make test` runs the testsets that hold this decision:

- `test/generators/superoperator.jl`: "Lindbladian channel rates are real and nonnegative by physical assumption"
- `test/generators/superoperator.jl`: "Lindbladian ingest rejects correlated channels and non-Hermitian drives"
- `test/generators/superoperator.jl`: "superoperator harmonics apply finite-level completeness"
- `test/generators/superoperator.jl`: "periodic Lindbladians lower through the public superoperator seam"
- `test/completion/state.jl`: "physical and lowered open-system constructors agree"
- `test/completion/gram.jl`: "automatic frame follows the canonical channel order"
- `test/completion/gram.jl`: "completed Gram accessors return independent containers"

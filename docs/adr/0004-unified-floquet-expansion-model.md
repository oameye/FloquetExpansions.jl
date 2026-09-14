# One result model and generic accessors

Hamiltonian and Liouvillian dynamics share `FloquetExpansion` and `floquet_expansion`. Consumers use `effective_generator`, `effective_component`, `micromotion`, `hamiltonian`, and `hamiltonian_component`; no separate Liouvillian result type is introduced.

Positive completion preserves this one-result model. `positive_completion(vv, algorithm)` returns another `FloquetExpansion` whose completion state records the selected finite positive continuation. The retained Van Vleck coefficients, gauge, and micromotion remain the same perturbative data; only the finite effective-generator realization changes.

A completed `FloquetExpansion` exposes representation-independent physical accessors such as `effective_generator`, `hamiltonian`, and `channels`, together with Kossakowski coordinates in its stored `DissipativeFrame`. Algorithm-specific intermediate data are available through a single expert `factorization` accessor rather than through a separate wrapper result type.

## Gate

`make test` runs the testsets that hold this decision:

- `test/completion_state.jl`: "positive-completion algorithms establish the public dispatch boundary"
- `test/completion_storage.jl`: "completed representation owns frame and cached retained data"
- `test/gram_completion.jl`: "completed Gram accessors return independent containers"
- `test/spectral_completion.jl`: "spectral factorization accessor returns independent containers"

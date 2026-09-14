# Raw finite-order effective generator

The package returns the direct finite-order van Vleck result for both Hamiltonian and Liouvillian inputs: an effective generator and micromotion. It does not reinterpret or replace a truncation when it no longer has Lindblad form; physical completion is a separate concern (see ADR 0006).

## Gate

`make test` runs the testsets that hold this decision:

- `test/engine.jl`: "the spec's closed forms, orders 0 to 2"
- `test/engine.jl`: "truncation follows the spec, X^[N] = sum_{k<N}"
- `test/completion_state.jl`: "raw Floquet expansions carry uncompleted state"
- `test/residual_scaling.jl`: "the truncated factorization solves eq:defining0 to O(wd^-(N-1))"

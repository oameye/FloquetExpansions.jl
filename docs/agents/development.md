# Development Workflow

Use this guide for source, test, documentation, or dependency changes.

## Before editing

1. Read the root `AGENTS.md` and check `git status --short --branch`. Preserve work that predates the task.
2. Read `CONTEXT.md` and the ADRs relevant to the behavior or representation being changed.
3. Identify the public API seam first. Tests should describe the behavior a package user can observe, so improve that seam when the existing API cannot express the required assertion.

The v0.0.1 policy permits breaking changes. Update every in-repository caller, test, docstring, example, and API page affected by a changed name, signature, or representation. Record a durable domain or design decision in `CONTEXT.md` or a new ADR when the change establishes project vocabulary or a module boundary.

## Package layout

- `src/FloquetExpansions.jl` wires modules and exports the public API.
- `src/periodic_operator.jl` owns the Fourier `PeriodicGenerator` boundary and harmonic calculus.
- `src/liouvillian.jl` owns the symbolic left/right Liouvillian representation, channels, composition, and Liouvillian harmonic lowering.
- `src/engine.jl` owns the generic Van Vleck recursion and `FloquetExpansion` accessors.
- `src/quasienergy.jl` owns symbolic Sambe-space `QuasienergyOperator` blocks.
- `test/*.jl` contains behavior tests; `test/helpers/` contains shared test support; `test/quality/` contains package-wide quality checks.
- `docs/src/` contains user-facing documentation, while `docs/adr/` contains design decisions and `docs/agents/` contains agent operating guidance.

Keep the shared Fourier and expansion seams generic across Hamiltonian and Liouvillian components. Keep sparse Liouvillian storage and numerical vectorization behind their respective interfaces.

## Test and quality loop

Use the smallest check that answers the current question, then run the full relevant gates before handoff:

- `make format` applies the repository's JuliaFormatter configuration.
- `make test` runs the package test suite and quality checks discovered by `test/runtests.jl`.
- `make jet` runs the separate JET audit.
- `make docs` checks docstrings and builds the documentation.
- `make bench` measures benchmark changes; report precompilation, TTFX, and steady-state runtime separately.

Run `make docs` from the repository root with Julia's normal depot and project
environment. Preserve the user's existing `JULIA_DEPOT_PATH` and `JULIA_PROJECT`
so Julia can reuse the established package downloads and precompiled caches. A
first run after a Julia, dependency, or source change may still precompile.

When the execution environment cannot write the normal depot, request writable
access or hand the command back to the user. Keep the validation command as
`make docs`; use a temporary depot only when the user explicitly requests an
isolated environment, because it discards cache reuse and can force a full
recompilation.

For new behavior, add a public-API regression test and a doctest or API documentation update when the user-facing interface changes. Use `@inferred` for type-stability contracts and `@allocations` for hot-path allocation contracts. Finish when the relevant checks pass and the diff accounts for the changed behavior, tests, and docs.

Generated documentation, local manifests, test logs, and benchmark output are ignored by the repository; leave them untracked and out of handoff summaries.

# Development Workflow

Use this guide for source, test, documentation, or dependency changes. It is the authority for the test layout, testing patterns, gate loop, and generated-documentation workflow. `STANDARDS.md` routes here.

## Before editing

1. Read the root `AGENTS.md` and check the working tree before changing it. Preserve work that predates the task.
2. Read `CONTEXT.md` and the ADRs relevant to the behavior or representation being changed.
3. Identify the public API seam first. Behavior tests should exercise what a package user can observe; improve that seam when the existing API cannot express the required assertion.

The v0.0.1 policy permits breaking changes. Update every in-repository caller, test, docstring, example, and API page affected by a changed name, signature, or representation. Record a durable domain or design decision in `CONTEXT.md` or an ADR when the change establishes project vocabulary or a module/representation boundary.

## Layout

[`architecture.md`](architecture.md) is the authority for what each file under `src/` owns. Around it:

- `test/*.jl` holds behavior tests, `test/helpers/` shared fixtures, and `test/quality/` package-wide checks.
- `docs/src/` holds user-facing documentation, `docs/adr/` design decisions, and `docs/agents/` repository guidance.
- `examples/*.jl` is the source for the tracked Literate pages under `docs/src/examples/`.
- `code_ratchet/` holds the CodeRatchet environment, hand-written rulings, and generated baseline files. Its configured measured scope is `src/`; other tracked Julia-code areas are explicitly classified as unmeasured in `rulings.toml`.

## Test structure

`test/runtests.jl` discovers tests with ParallelTestRunner. In the normal unfiltered `Pkg.test()` path it deliberately removes `quality/JET`, because JET is a separate acceptance gate with its own environment and workflow. Adding an ordinary `test/**/*.jl` file is otherwise enough to register it; there is no central test list to update.

Passing positional test arguments filters the discovered suite. Keep JET invocation explicit through `make jet` rather than relying on the default test command.

The qualified expert API is marked with `@public` in `src/FloquetExpansions.jl` rather than exported. `test/runtests.jl` explicitly imports those names into each isolated test worker. When that expert seam changes, update the module declaration and the test-worker import list together.

## Testing patterns

- **Test user-visible behavior through the public API.** Behavior and regression tests for package features should exercise exported or intentionally qualified `@public` interfaces. Do not reach through private fields or helpers merely to make such a test convenient.
- **Use internal tests for internal invariants.** Focused algebra or algorithm tests may qualify private helpers when directly validating a load-bearing invariant that is not usefully observable at the public seam, as in matrix-series or Gram-recursion tests. Such a test does not make that helper public API.
- Use `@inferred` when a stable return-type contract is part of the behavior.
- Keep compiler-sensitive core workloads under `JET.@test_opt` when optimizer cleanliness is an acceptance property. The current completion workloads in `test/quality/JET.jl` are required gates, not optional diagnostics.
- For runtime-sensitive changes, measure the benchmark workload that exercises the path. Add an allocation assertion only when a small, stable operation is genuinely expected to have a fixed allocation contract; the repository does not currently impose a package-wide zero-allocation gate.
- For new user-facing behavior, add a public-API regression test and update a doctest or manual/API documentation where that interface is documented.

## The gate loop

Use the smallest check that answers the current question, then run the relevant acceptance gates before handoff. **`make test` and `make jet` are separate gates and neither covers the other.**

| Command | Runs | CI |
| --- | --- | --- |
| `make format` | JuliaFormatter in-place over the repository | `Format.yml` checks formatting |
| `make test` | the default ParallelTestRunner suite, including Aqua, CheckConcreteStructs, ExplicitImports and doctests, but excluding `quality/JET` | `Tests.yml` |
| `make jet` | `test/quality/JET.jl`: package JET plus explicit optimizer-stability workloads | `JET.yml` |
| `make ratchet` | configured CodeRatchet metrics; locally includes the JET metric | `Ratchet.yml`, with duplicate JET omitted there because `JET.yml` is stronger |
| `make docs` | the Documenter build with `checkdocs=:exports`; doctests are disabled here because the test suite owns them | `Documentation.yml` |
| `make bench` | the benchmark suite | `Benchmarks.yaml` |

`make all` is `setup format test docs`; it does not run JET, the ratchet, benchmarks, or the CI-only spell check.

Spelling is checked on pull requests by `SpellCheck.yml`, configured by `.typos.toml`.

## Building the docs

Run `make docs` from the repository root with Julia's normal depot and project environment so existing downloads and precompile caches can be reused. A first run after a Julia, dependency, or source change may still precompile.

`docs/make.jl` sets `doctest=false` because doctests run in `test/quality/Documenter.jl`. It uses `checkdocs=:exports` for exported-doc inclusion.

## Generated files

`examples/*.jl` is the source of truth for `docs/src/examples/*.md`. Under documentation CI, `docs/make.jl` includes `docs/make_md_examples.jl`, which runs Literate before `makedocs`; a normal local `make docs` without the CI condition builds the tracked Markdown as it stands and does not regenerate it. Do not edit the generated Markdown as the source of a documentation change.

`docs/build/`, `docs/site/`, `Manifest.toml`, `test-run.log`, and `benchmark/benchmarks_output.json` are gitignored. CodeRatchet refresh/triage output under `code_ratchet/_refresh/` and `code_ratchet/_triage/` is also generated and ignored. Leave generated local artifacts out of commits and handoff summaries.

## Finishing

Finish when the gates for the area you changed pass and the diff accounts for the changed behavior, its tests, and its docs. Run `make ratchet` for `src/` changes and `make bench` when the change is performance-sensitive.

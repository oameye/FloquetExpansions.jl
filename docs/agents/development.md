# Development Workflow

Use this guide for source, test, documentation, or dependency changes. Authority for the test layout, the testing patterns, and the gate loop. `STANDARDS.md` routes here.

## Before editing

1. Read the root `AGENTS.md` and check `git status --short --branch`. Preserve work that predates the task.
2. Read `CONTEXT.md` and the ADRs relevant to the behavior or representation being changed.
3. Identify the public API seam first. Tests describe the behavior a package user can observe, so improve that seam when the existing API cannot express the required assertion.

The v0.0.1 policy permits breaking changes. Update every in-repository caller, test, docstring, example, and API page affected by a changed name, signature, or representation. Record a durable domain or design decision in `CONTEXT.md` or a new ADR when the change establishes project vocabulary or a module boundary.

## Layout

[`architecture.md`](architecture.md) is the authority for what each file under `src/` owns; read its Module ownership table rather than a second list here. Around it:

- `test/*.jl` holds behavior tests, `test/helpers/` shared fixtures, `test/quality/` the package-wide checks.
- `docs/src/` holds user-facing documentation, `docs/adr/` design decisions, `docs/agents/` agent operating guidance.

## Test structure

`test/runtests.jl` discovers every `test/**/*.jl` and runs it with ParallelTestRunner. Adding a file is all it takes to register a test; there is no list to update.

Passing a positional argument selects tests by prefix, so `Pkg.test(; test_args=["quality"])` runs the quality checks alone.

## Testing patterns

- **Test through the public API.** Exercise the package as a user would, through exported or documented interfaces. Where behavior cannot be reached that way, improve the public seam rather than adding a test-only escape hatch. This one is held by review; no check enforces it.
- Assert type-stability contracts with `@inferred`.
- Assert hot-path allocation contracts with `@allocations`. See [`performance.md`](performance.md) for which paths those are.
- For new behavior, add a public-API regression test, plus a doctest or API documentation update when the user-facing interface changes.

## The gate loop

Use the smallest check that answers the current question, then run the gates below before handoff. **`make test` and `make jet` are separate gates and neither covers the other.**

| Command | Runs | CI |
| --- | --- | --- |
| `make format` | JuliaFormatter over the repository | `Format.yml` |
| `make test` | every test except `quality/JET`, including Aqua, CheckConcreteStructs, ExplicitImports and the doctests | `Tests.yml` |
| `make jet` | `test/quality/JET.jl` alone | `JET.yml` |
| `make docs` | docstring checks and the documentation build | `Documentation.yml` |
| `make bench` | the benchmark suite | `Benchmarks.yaml` |

`test/runtests.jl` deletes `quality/JET` whenever no positional argument is given, which is exactly the case `make test` produces. So a green `make test` says nothing about JET, and `make jet` is the only local command that does. `make all` is `setup format test docs`, so it does not run JET either despite the name.

Spelling is checked in CI only, by `SpellCheck.yml` over the whole tree, configured by `.typos.toml`.

## Building the docs

Run `make docs` from the repository root with Julia's normal depot and project environment. Preserve the existing `JULIA_DEPOT_PATH` and `JULIA_PROJECT` so Julia reuses established downloads and precompiled caches. A first run after a Julia, dependency, or source change may still precompile.

Where the environment cannot write the normal depot, request writable access or hand the command back to the user. Keep the validation command as `make docs`. Use a temporary depot only on explicit request, because it discards cache reuse and can force a full recompilation.

## Generated files

`docs/make_md_examples.jl` runs Literate over `examples/*.jl` and writes `docs/src/examples/*.md` on every docs build. **Those `.md` files are tracked in git but regenerated, so an edit to one is overwritten.** Edit the `examples/*.jl` source instead.

`docs/build/`, `Manifest.toml`, `test-run.log` and `benchmark/benchmarks_output.json` are gitignored. Leave them untracked and out of handoff summaries.

## Finishing

Finish when the gates for the area you changed pass, and the diff accounts for the changed behavior, its tests, and its docs.

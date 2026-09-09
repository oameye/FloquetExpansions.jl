# FloquetExpansions.jl Standards Map

This file is a **map, not a rule**. It carries no standards of its own: it tells you which file governs the thing you are about to change, and which check fails if you get it wrong.

Read the row for your subject, open the file it names, and follow that file. Where no row covers what you are doing, the standard does not exist yet. Say so rather than inventing one.

This map's structure is adapted from [PortfolioOptimisers.jl](https://github.com/dcelisgarza/PortfolioOptimisers.jl)
by Daniel Celis Garza (MIT), described in [this Discourse post](https://discourse.julialang.org/t/idiomatic-julia-code-in-ai-generated-code/139183/4).
The map-not-rule framing, the Rule/Scope/Gate/Authority vocabulary, the precedence ladder and the three-table layout are his.

## Vocabulary

- **Rule**: one normative statement in a standards file.
- **Scope**: the files or behavior a rule governs.
- **Gate**: an automated check that fails on a breach. A rule with no gate is held by review.
- **Authority**: the one file that owns a rule's text. Other files link to it instead of maintaining a second copy.

## Precedence

Where two files disagree, the higher entry wins. Report the disagreement rather than silently picking a side: a contradiction between standards files is itself a defect.

1. **`docs/adr/`**: a decision that reached `main` outranks every other file on the point it settles. An ADR describing superseded behaviour is correct history, not a bug.
2. **`CONTEXT.md`**: the domain glossary. It fixes the project vocabulary.
3. **`AGENTS.md`**: repository-wide working agreements.
4. **`docs/agents/*.md`**: per-scope guides.
5. **`README.md`** and `docs/src/`: user-facing documentation, which describes behaviour rather than governing it.

## What am I about to touch?

| Subject | Authority | Gate |
| --- | --- | --- |
| A struct field's storage type | [`docs/agents/style.md`](docs/agents/style.md) | `test/quality/CheckConcreteStructs.jl` via `make test` |
| An import or qualified dependency access | [`docs/agents/style.md`](docs/agents/style.md) | `test/quality/ExplicitImports.jl` via `make test` |
| Formatting | `.JuliaFormatter.toml` | `.github/workflows/Format.yml` |
| A function signature's type constraints | [`docs/agents/style.md`](docs/agents/style.md) | none, unenforced |
| Keyword syntax or forwarding form | [`docs/agents/style.md`](docs/agents/style.md) | none, unenforced |
| A comment, or a docstring on an internal | [`docs/agents/style.md`](docs/agents/style.md) | none, unenforced |
| A name for a private internal | [`docs/agents/style.md`](docs/agents/style.md) | none, unenforced |
| A name for a domain concept | [`CONTEXT.md`](CONTEXT.md) | none, unenforced |
| An `export` line | `src/FloquetExpansions.jl` | `test/quality/Aqua.jl` for undefined exports only |
| A qualified expert `@public` name | `src/FloquetExpansions.jl` | none, unenforced |
| Type inference on a package or core path | [`docs/agents/performance.md`](docs/agents/performance.md) | `make jet`, `.github/workflows/JET.yml`, plus targeted `@inferred` / `JET.@test_opt` tests where the contract is load-bearing |
| Optimizer stability of the explicit completion core | [`docs/agents/performance.md`](docs/agents/performance.md) | `test/quality/JET.jl` contains required `JET.@test_opt` workloads |
| Runtime performance regression | [`docs/agents/performance.md`](docs/agents/performance.md) | `.github/workflows/Benchmarks.yaml`: alert above 130%, fail above 170% of baseline |
| Allocation behaviour | [`docs/agents/performance.md`](docs/agents/performance.md) | none repository-wide; measure explicitly and add a targeted test only where a stable allocation contract exists |
| How a performance number is reported | [`docs/agents/performance.md`](docs/agents/performance.md) | none, unenforced |
| Adding or changing a test file | [`docs/agents/development.md`](docs/agents/development.md) | `test/runtests.jl` auto-discovers tests; the default suite deliberately excludes `quality/JET` |
| Reaching past the public API in a behavior test | [`docs/agents/development.md`](docs/agents/development.md) | none, unenforced |
| A docstring on an exported or qualified expert name | [`docs/agents/development.md`](docs/agents/development.md) | none, unenforced |
| A `jldoctest` block | [`docs/agents/development.md`](docs/agents/development.md) | `test/quality/Documenter.jl` via `make test` |
| Inclusion of exported docstrings in the manual | `docs/make.jl` (`checkdocs=:exports`) | `make docs`, `.github/workflows/Documentation.yml` |
| A module boundary or symbolic representation | [`docs/agents/architecture.md`](docs/agents/architecture.md), and the ADR it cites | none, unenforced unless a behavior test covers the consequence |
| Positive-completion result and representation semantics | [`docs/adr/0009-cp-preserving-floquet-completion.md`](docs/adr/0009-cp-preserving-floquet-completion.md) | `make test`, especially the completion-state/storage/validation tests |
| The dissipative quasienergy convention | [`docs/adr/0008-dissipative-quasienergy-operator.md`](docs/adr/0008-dissipative-quasienergy-operator.md) | `make test` |
| A dependency or a `[compat]` bound | `Project.toml` | `test/quality/Aqua.jl` via `make test` |
| A file under `docs/src/examples/` | [`docs/agents/development.md`](docs/agents/development.md) § Generated files | none; documentation CI regenerates it from `examples/*.jl` |
| Prose anywhere in the repository | `.typos.toml` | `.github/workflows/SpellCheck.yml` on pull requests |
| A breaking change's blast radius | [`AGENTS.md`](AGENTS.md) § Development stage | none, unenforced |
| Committing or pushing | [`AGENTS.md`](AGENTS.md) § Git policy | none, unenforced |
| A decision worth recording | [`docs/adr/`](docs/adr/) | none, unenforced |

## The standards files

| File | Owns | Scope |
| --- | --- | --- |
| `STANDARDS.md` | this routing map | the repository |
| `CONTEXT.md` | project vocabulary | the repository |
| `AGENTS.md` | orientation, development stage, git policy | the repository |
| `docs/adr/` | durable design decisions and their reasoning | named per ADR |
| `docs/agents/development.md` | repository workflow, test layout and patterns, gate loop, generated documentation | repository changes; especially `test/` and `docs/` |
| `docs/agents/architecture.md` | data flow, source-file ownership, public seams, representation rules | `src/` |
| `docs/agents/style.md` | signatures, fields, imports, formatting, comments, names | `src/` |
| `docs/agents/performance.md` | inference, optimizer stability, allocations, hot-path shape, performance terminology | `src/`, `benchmark/` |
| `docs/agents/domain.md` | how project vocabulary and ADRs are consumed | the repository |
| `docs/agents/issue-tracker.md` | issue and spec conventions | GitHub Issues |
| `docs/agents/triage-labels.md` | the label set | GitHub Issues |

## The gates

Every gate below is described by what the current repository actually executes.

| Gate | Enforces | How to run |
| --- | --- | --- |
| `test/quality/CheckConcreteStructs.jl` | `all_concrete(FloquetExpansions)`; no abstract/`Any` storage fields in package structs | `make test`, or select the quality test |
| `test/quality/ExplicitImports.jl` | no implicit imports, no stale explicit imports, owner-correct imports/qualified accesses, no self-qualified accesses | `make test`, or select the quality test |
| `test/quality/Aqua.jl` | Aqua's package checks through `Aqua.test_all(...; project_extras=false)` | `make test`, or select the quality test |
| `test/quality/Documenter.jl` | package doctests | `make test`, or select the quality test |
| `test/quality/JET.jl` | package-wide JET reports plus explicit `JET.@test_opt` workloads for Gram, recursive Gram, Spectral, and completed accessors | **`make jet` only.** See below |
| `test/runtests.jl` | test discovery and execution through ParallelTestRunner; the unfiltered default removes `quality/JET` | `make test` |
| `.github/workflows/Tests.yml` | the default `Pkg.test()` suite on supported CI Julia | CI when its path filter matches |
| `.github/workflows/JET.yml` | `test/quality/JET.jl` | CI when its path filter matches |
| `.github/workflows/Format.yml` | JuliaFormatter with repository configuration over the tree | `make format`, then `jlfmt --check --verbose .` for a non-mutating check |
| `.github/workflows/Benchmarks.yaml` | runtime benchmark non-regression against stored benchmark history | `make bench` |
| `.github/workflows/SpellCheck.yml` | spelling, configured by `.typos.toml` | CI only |
| `.github/workflows/Documentation.yml` | documentation build and exported-doc inclusion; CI also regenerates Literate examples before `makedocs` | `make docs` for the local build |

**`make test` does not run JET.** `test/runtests.jl` removes `quality/JET` from the unfiltered suite used by `Pkg.test()`, and therefore by `make test`. A green `make test` says nothing about the dedicated optimizer/inference acceptance workloads. Run `make jet` separately. In CI they are separate workflows, `Tests.yml` and `JET.yml`.

`make all` is `setup format test docs`; despite its name it does not include JET, benchmarks, or the CI-only spell check.

Several CI workflows are path-filtered. If a relevant workflow is skipped because the pull request did not touch its configured paths, that skip is not evidence that another local gate ran in its place.

A row above reading **none, unenforced** names a rule that no current gate checks. That is a known state rather than an implied guarantee.

## Changing a standard

Amending a rule is a deliberate change, not a passing edit made while fixing something else.

1. Find the **authority** for the rule in the tables above. Change the text there, and nowhere else.
2. Where another file restates the rule, replace the copy with a link. Two live copies drift.
3. Update the row in this file when the authority, scope, or gate changed.
4. Where a load-bearing rule has no gate, consider adding one that tests the actual contract rather than a proxy.
5. Where the change reverses a durable design decision, update or supersede the relevant ADR.

Audit the map both ways: current repository behavior must be routed by the map where a standard exists, and every statement in the map must remain true of the repository.

## Before you finish a change

1. `make format` has been applied and the formatter check is clean.
2. `make test` passes for source, test, or dependency changes.
3. `make jet` passes when `src/`, inference-sensitive dependencies, or compiler-sensitive tests changed.
4. `make docs` passes when docstrings or `docs/` changed.
5. `make bench` is run for performance-sensitive source changes.
6. The diff accounts for the changed behaviour, its tests, and its documentation.

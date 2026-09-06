# FloquetExpansions.jl Standards Map

This file is a **map, not a rule**. It carries no standards of its own: it tells you which file governs the thing you are about to change, and which check fails if you get it wrong.

Read the row for your subject, open the file it names, and follow that file. Where no row covers what you are doing, the standard does not exist yet. Say so rather than inventing one.

This map's structure is adapted from [PortfolioOptimisers.jl](https://github.com/dcelisgarza/PortfolioOptimisers.jl)
by Daniel Celis Garza (MIT), described in [this Discourse post](https://discourse.julialang.org/t/idiomatic-julia-code-in-ai-generated-code/139183/4).
The map-not-rule framing, the Rule/Scope/Gate/Authority vocabulary, the precedence ladder and the three-table layout are his.

## Vocabulary

Four words, used precisely here and in the `improve-codebase-maintainability` skill that audits this map.

- **Rule**: one normative statement in a standards file.
- **Scope**: the files a rule governs. A scope is real only if you can run it and count.
- **Gate**: the automated check that fails when the rule breaks. A rule with no gate holds only while every contributor remembers it.
- **Authority**: the one file that owns a rule's text. Every other mention links here instead of restating it, because two copies drift.

## Precedence

Where two files disagree, the higher entry wins. Report the disagreement rather than silently picking a side: a contradiction between standards files is itself a defect.

1. **`docs/adr/`**: a decision that reached `main` outranks every other file on the point it settles. An ADR describing superseded behaviour is correct history, not a bug.
2. **`CONTEXT.md`**: the domain glossary. It fixes the words; nothing else may rename a concept.
3. **`AGENTS.md`**: the working agreements for this checkout.
4. **`docs/agents/*.md`**: the per-scope guides.
5. **`README.md`** and `docs/src/`: user-facing documentation, which describes behaviour rather than governing it.

## What am I about to touch?

| Subject | Authority | Gate |
| --- | --- | --- |
| A struct field's type | [`docs/agents/style.md`](docs/agents/style.md) | `test/quality/CheckConcreteStructs.jl` |
| An import line | [`docs/agents/style.md`](docs/agents/style.md) | `test/quality/ExplicitImports.jl` |
| Formatting | `.JuliaFormatter.toml` | `.github/workflows/Format.yml` |
| A function signature's types | [`docs/agents/style.md`](docs/agents/style.md) | none, unenforced |
| A keyword argument's form | [`docs/agents/style.md`](docs/agents/style.md) | none, unenforced |
| A comment, or a docstring on an internal | [`docs/agents/style.md`](docs/agents/style.md) | none, unenforced |
| A name for a private internal | [`docs/agents/style.md`](docs/agents/style.md) | none, unenforced |
| A name for a domain concept | [`CONTEXT.md`](CONTEXT.md) | none, unenforced |
| Type inference on a new path | [`docs/agents/performance.md`](docs/agents/performance.md) | `make jet`, `.github/workflows/JET.yml`, plus `@inferred` in the behaviour's test |
| Allocations on a core path | [`docs/agents/performance.md`](docs/agents/performance.md) | `@allocations` in the test; `.github/workflows/Benchmarks.yaml` fails above 170% |
| How a performance number is reported | [`docs/agents/performance.md`](docs/agents/performance.md) | none, unenforced |
| Adding or changing a test file | [`docs/agents/development.md`](docs/agents/development.md) | `test/runtests.jl` auto-discovery |
| Reaching past the public API in a test | [`docs/agents/development.md`](docs/agents/development.md) | none, unenforced |
| A docstring on a public name | [`docs/agents/development.md`](docs/agents/development.md) | `test/quality/Documenter.jl` |
| A `jldoctest` block | [`docs/agents/development.md`](docs/agents/development.md) | `test/quality/Documenter.jl` |
| A module boundary or a representation | [`docs/agents/architecture.md`](docs/agents/architecture.md), and the ADR it cites | none, unenforced |
| An `export` line | `src/FloquetExpansions.jl` | `test/quality/Aqua.jl` (undefined exports only) |
| A dependency or a `[compat]` bound | `Project.toml` | `test/quality/Aqua.jl` |
| A file under `docs/src/examples/` | [`docs/agents/development.md`](docs/agents/development.md) § Generated files | none; `make docs` overwrites the file |
| Prose anywhere in the repository | `.typos.toml` | `.github/workflows/SpellCheck.yml` |
| A breaking change's blast radius | [`AGENTS.md`](AGENTS.md) § Development stage | none, unenforced |
| Committing or pushing | [`AGENTS.md`](AGENTS.md) § Git policy | none, unenforced |
| A decision worth recording | [`docs/adr/`](docs/adr/) | none, unenforced |

## The standards files

| File | Owns | Scope |
| --- | --- | --- |
| `STANDARDS.md` | this map | the repository |
| `CONTEXT.md` | the domain glossary, one or two sentences per term | the repository |
| `AGENTS.md` | orientation, development stage, git policy | the repository |
| `docs/adr/` | decisions and their reasoning, one file per decision | named per ADR |
| `docs/agents/development.md` | the workflow, the test layout, the testing patterns, the gate loop, generated files | `test/`, `docs/` |
| `docs/agents/architecture.md` | the data flow, module ownership, representation rules | `src/**/*.jl` |
| `docs/agents/style.md` | signatures, types, imports, formatting, comments, names | `src/**/*.jl` |
| `docs/agents/performance.md` | inference, allocations, hot-path shape, the three axes | `src/**/*.jl`, `benchmark/` |
| `docs/agents/domain.md` | how the skills consume `CONTEXT.md` and `docs/adr/` | the repository |
| `docs/agents/issue-tracker.md` | issue and spec conventions | GitHub Issues |
| `docs/agents/triage-labels.md` | the label set | GitHub Issues |

## The gates

Every gate below fails on a real breach.

| Gate | Enforces | How to run |
| --- | --- | --- |
| `test/quality/CheckConcreteStructs.jl` | every struct field is concretely typed | `make test`, or the file |
| `test/quality/ExplicitImports.jl` | no implicit import, no stale explicit import, no non-owner import or qualified access, no self-qualified access | `make test`, or the file |
| `test/quality/Aqua.jl` | method ambiguities, unbound type parameters, undefined exports, stale dependencies, `[compat]` completeness, type piracy, persistent tasks | `make test`, or the file |
| `test/quality/Documenter.jl` | every `jldoctest` still produces its printed output | `make test`, or the file |
| `test/quality/JET.jl` | JET reports nothing in `FloquetExpansions` | **`make jet` only.** See below |
| `test/runtests.jl` | every `test/**/*.jl` is discovered and run | `make test` |
| `.github/workflows/Format.yml` | `blue` style, indent 2, over the whole repository | `make format`, then `jlfmt --check .` |
| `.github/workflows/Benchmarks.yaml` | no benchmark above 170% of baseline; comments above 130% | `make bench` |
| `.github/workflows/SpellCheck.yml` | spelling, minus the `.typos.toml` allow-list | CI only |
| `.github/workflows/Documentation.yml` | the documentation builds | `make docs` |

**`make test` does not run JET.** `test/runtests.jl` deletes `quality/JET` from the suite whenever no positional argument is given, which is the case `Pkg.test()` and therefore `make test` produce. A green `make test` says nothing about JET, and `make jet` is the only local command that runs it. In CI they are separate workflows, `Tests.yml` and `JET.yml`.

A row in the previous table reading **none, unenforced** names a rule that no gate checks. That is a known state, not a hidden one: an unenforced rule holds by review and by memory. Roughly half the rows are in that state; count them in the table rather than trusting a number written here, because a written count goes stale where it stands.

The largest cluster of them is the representation rules in `architecture.md`, and those are the ones most expensive to get wrong.

## Changing a standard

Amending a rule is the maintainer's call. It is never a passing edit made while fixing something else.

1. Find the **authority** for the rule in the tables above. Change the text there, and nowhere else.
2. Where another file restates the rule, replace the copy with a link. Two live copies drift.
3. Update the row in this file when the authority, the scope, or the gate changed.
4. Where a load-bearing rule has no gate, consider adding one. Keep it quiet under ordinary work, and open it with a comment saying what drifted, when, and why the check earns its runtime.
5. Where the change reverses a previous decision, record it in `docs/adr/`.

The `improve-codebase-maintainability` skill audits this map both ways: the code against each rule, and each rule against the code.

## Before you finish a change

1. `make format` passes.
2. `make test` passes.
3. `make jet` passes, when you touched `src/`.
4. `make docs` passes, when you touched a docstring or `docs/`.
5. The diff accounts for the changed behaviour, its tests, and its docs.

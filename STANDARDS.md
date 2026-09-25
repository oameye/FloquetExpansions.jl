# FloquetExpansions.jl Standards Map

A map, not a rule. Every rule's text lives in the file named under **Authority**. This file adds
one thing only: which check fails when the rule breaks. Where no row covers what you are about to
change, the standard does not exist yet. Say so rather than inventing one.

## Precedence

Where two files disagree the higher entry wins. Report the disagreement rather than picking a
side: a contradiction between standards files is itself a defect.

1. `docs/adr/`, for the point each ADR settles. An ADR describing superseded behaviour is correct
   history, not a bug. Each ADR names its own gate, so no ADR is routed from here.
2. `CONTEXT.md`, the domain glossary.
3. `AGENTS.md`, repository-wide working agreements.
4. `docs/agents/*.md`, per-scope guides.
5. `README.md` and `docs/src/`, which describe behaviour rather than govern it.

## Gated rules

| Subject | Authority | Gate |
| --- | --- | --- |
| Struct field storage types | [`docs/agents/style.md`](docs/agents/style.md) | `test/quality/CheckConcreteStructs.jl`, via `make test` |
| Imports and qualified access | [`docs/agents/style.md`](docs/agents/style.md) | `test/quality/ExplicitImports.jl`, via `make test` |
| Public-API-only package tests | [`docs/agents/development.md`](docs/agents/development.md) | `test/quality/PublicAPI.jl`, via `make test` |
| Formatting | `.JuliaFormatter.toml` | `.github/workflows/Format.yml`; `make format` locally |
| A dependency or `[compat]` bound | `Project.toml` | `test/quality/Aqua.jl`, via `make test` |
| An `export` line | `src/FloquetExpansions.jl` | `test/quality/Aqua.jl`, undefined exports only |
| A `jldoctest` block | [`docs/agents/development.md`](docs/agents/development.md) | `test/quality/Documenter.jl`, via `make test` |
| Exported docstrings reaching the manual | `docs/make.jl` (`checkdocs=:exports`) | `.github/workflows/Documentation.yml`; `make docs` |
| Type inference and optimizer stability | [`docs/agents/performance.md`](docs/agents/performance.md) | `make jet` and `.github/workflows/JET.yml`. `make test` excludes it |
| Runtime performance | [`docs/agents/performance.md`](docs/agents/performance.md) | `.github/workflows/Benchmarks.yaml`: alert at 130%, fail at 170%; `make bench` |
| Prose spelling | `.typos.toml` | `.github/workflows/SpellCheck.yml`, on pull requests |

## Held by review alone

No gate checks the following. That is a known state, not an implied guarantee.

| Subject | Authority |
| --- | --- |
| Signature type constraints, keyword and forwarding form, comments, internal docstrings, private names | [`docs/agents/style.md`](docs/agents/style.md) |
| Allocation behaviour, and how a performance number is reported | [`docs/agents/performance.md`](docs/agents/performance.md) |
| Semantic private-field access that cannot be recognized statically, and docstrings on qualified `@public` names | [`docs/agents/development.md`](docs/agents/development.md) |
| Module boundaries, symbolic representation, and the qualified `@public` seam | [`docs/agents/architecture.md`](docs/agents/architecture.md) |
| A name for a domain concept | [`CONTEXT.md`](CONTEXT.md) |
| A breaking change's blast radius, and git policy | [`AGENTS.md`](AGENTS.md) |
| Whether a decision gets recorded | [`docs/adr/`](docs/adr/) |

## Maintaining this file

Change a rule's text in its authority and nowhere else. Update a row here only when the authority
or the gate changed, not when the rule's wording did.

Several CI workflows are path-filtered. A workflow skipped because the pull request did not touch
its configured paths is not evidence that another gate ran in its place.

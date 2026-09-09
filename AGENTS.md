# FloquetExpansions.jl

A Julia package for symbolic high-frequency expansions of periodically driven quantum systems. The current public implementation is the Van Vleck expansion for Hamiltonian and Liouvillian generators.

## Orientation

Four files carry the rules, and each has one job.

- **[`STANDARDS.md`](STANDARDS.md) is the map.** It routes a subject to the file that owns the rule and to the check that fails when the rule breaks. Start there when you do not know which file governs what you are about to change. It carries no rules of its own.
- **[`CONTEXT.md`](CONTEXT.md) is the glossary.** It fixes the domain language and nothing else may rename a concept. Read it before asking a question. Introducing or renaming a term means updating it in the same change.
- **[`docs/adr/`](docs/adr/) records decisions** and the reasoning behind them, one file per decision. A decision that reached `main` outranks every other file on the point it settles.
- **[`docs/agents/`](docs/agents/) holds the per-scope guides.** `development.md` for the workflow and the gates, `architecture.md` for the module map, `style.md` and `performance.md` for how source is written, `domain.md`, `issue-tracker.md` and `triage-labels.md` for the rest.

One rule, one authority. Where two of these disagree, report it rather than picking a side: the contradiction is itself a defect.

## Development stage

The package is at v0.0.1 and its public API and internal representations are experimental. Breaking changes are welcome when they improve correctness, clarity, or the domain model. When making one, update all in-repository callers, tests, documentation, `CONTEXT.md`, and relevant ADRs in the same change. A compatibility alias or migration shim needs an explicit reason.

## Git policy

**Commits are the maintainer's to make.** Neither Claude nor any subagent runs `git commit`, `git push`, or any git command that modifies history, except when explicitly instructed otherwise.

## Working

For a source, test, documentation, or dependency change, read [`docs/agents/development.md`](docs/agents/development.md). It owns the gate loop, and the one thing worth knowing before reading it is that **`make test` and `make jet` are separate gates and neither covers the other.** Run `make help` for the full target list.

For module boundaries or representation changes, also read [`docs/agents/architecture.md`](docs/agents/architecture.md) and the relevant ADRs.

Use the `julia-mcp` server (`julia_eval`, `julia_list_sessions`, `julia_restart`) for quick questions: checking a type, evaluating an expression, verifying a method signature. It holds a persistent session, so prefer it over a full test run when you only need an answer.

## Issues

Issues and specs live in GitHub Issues for `oameye/FloquetExpansions.jl`; use the `gh` CLI. See [`docs/agents/issue-tracker.md`](docs/agents/issue-tracker.md), and [`docs/agents/triage-labels.md`](docs/agents/triage-labels.md) for the labels.

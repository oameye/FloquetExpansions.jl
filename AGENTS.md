# FloquetExpansions.jl

A Julia package for symbolic high-frequency expansions of periodically driven quantum systems. The current public expansion method is Van Vleck for Hamiltonian and Liouvillian generators. Liouvillian expansions can additionally be passed through an explicit positive-completion stage using `Gram()` or the restricted `Spectral()` realization.

## Orientation

Four entry points organize the repository guidance.

- **[`STANDARDS.md`](STANDARDS.md) is the map.** It routes a subject to the file that owns the rule and to the check that enforces it. Start there when you do not know which file governs what you are about to change. It carries no rules of its own.
- **[`CONTEXT.md`](CONTEXT.md) is the glossary.** It fixes the domain language. Introducing or renaming a domain concept means updating it in the same change.
- **[`docs/adr/`](docs/adr/) records decisions** and the reasoning behind them, one file per decision. A decision that reached `main` outranks the working guides on the point it settles.
- **[`docs/agents/`](docs/agents/) holds the per-scope guides.** `development.md` owns workflow and gates, `architecture.md` module and representation boundaries, and `style.md` / `performance.md` source-writing rules. `domain.md`, `issue-tracker.md`, and `triage-labels.md` cover the remaining repository workflows.

One rule, one authority. Where two authorities disagree, report the contradiction rather than silently choosing one.

## Development stage

The package is at v0.0.1 and its public API and internal representations are experimental. Breaking changes are welcome when they improve correctness, clarity, or the domain model. When making one, update all in-repository callers, tests, documentation, `CONTEXT.md`, and relevant ADRs in the same change. A compatibility alias or migration shim needs an explicit reason.

## Git policy

**Commits are the maintainer's to make.** Neither an agent nor a subagent runs `git commit`, `git push`, or another history-modifying git command unless explicitly instructed otherwise.

## Working

For a source, test, documentation, or dependency change, read [`docs/agents/development.md`](docs/agents/development.md). It owns the gate loop. The important distinction is that **`make test` and `make jet` are separate gates and neither covers the other.** Run `make help` for the repository's local targets.

For module boundaries, representation changes, completion semantics, or quasienergy conventions, also read [`docs/agents/architecture.md`](docs/agents/architecture.md) and the relevant ADRs.

## Issues

Issues and specs live in GitHub Issues for `oameye/FloquetExpansions.jl`. See [`docs/agents/issue-tracker.md`](docs/agents/issue-tracker.md) for the repository's issue workflow and [`docs/agents/triage-labels.md`](docs/agents/triage-labels.md) for labels.

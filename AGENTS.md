# FloquetExpansions.jl

## What is this?

A Julia package for symbolic high-frequency expansions of periodically driven quantum systems. The current public implementation is the Van Vleck expansion for Hamiltonian and Liouvillian generators.

## Development stage

The package is at v0.0.1 and its public API and internal representations are experimental. Breaking changes are welcome when they improve correctness, clarity, or the domain model. When making one, update all in-repository callers, tests, documentation, `CONTEXT.md`, and relevant ADRs in the same change. Compatibility aliases and migration shims need an explicit reason.

## Git policy

**Never commit or push.** Neither Claude nor any subagent may run `git commit`, `git push`, or any git command that modifies history. All commits are made by the user, except when explicitly instructed otherwise.

## Development workflow

For a source, test, documentation, or dependency change, read [`docs/agents/development.md`](docs/agents/development.md). For module boundaries or representation changes, also read [`docs/agents/architecture.md`](docs/agents/architecture.md) and the relevant ADRs.

Common tasks are defined by the Makefile:

```sh
make test       # run all tests (Pkg.test)
make format     # format with JuliaFormatter
make docs       # build documentation
make servedocs  # serve docs locally with LiveServer
make bench      # run benchmarks
make all        # setup + format + test + docs
```

### Quick debugging with Julia MCP

Use the `julia-mcp` MCP server (tools: `julia_eval`, `julia_list_sessions`, `julia_restart`) for quick debugging and testing small snippets — e.g., checking a type, evaluating an expression, or verifying a method signature. Prefer this over spinning up a full test run when you just need a quick answer.

### Test structure

Tests for each module are `test/*.jl`; shared fixtures live in `test/helpers/`, and static quality checks live in `test/quality/`. `test/runtests.jl` discovers and runs them with ParallelTestRunner.

### Testing patterns

- **Golden rule: tests use only the public API.** Exercise the package as a user
  would through exported or documented interfaces; do not call private methods,
  reach into internal submodules or fields, or assert implementation details.
  A good test exercises a feature through the public API and checks it against
  an expected result. If behavior cannot be tested through the public API,
  improve the public seam instead of adding a test-only escape hatch.
- `@inferred` for type stability checks
- `@allocations` for zero-allocation verification on hot paths

### Quality gates

Before merging any PR:
1. `make test` passes (all quality + unit + integration tests)
2. JET reports zero issues
3. No `Any`-typed fields in any struct (CheckConcreteStructs)
4. All imports explicit (ExplicitImports)

The repository uses JuliaFormatter, configured by `.JuliaFormatter.toml` and
invoked through `make format`.

## Coding rules

### Function signatures

- **Use the most restrictive signature type possible.** Tight type declarations let JET catch errors. When prototyping it's fine to start loose, but committed code should have specific types.
- **Explicit `;` for keyword arguments.** Always use an explicit semicolon:
  ```julia
  # Good
  FockSpace(; name = :a)
  # Bad
  FockSpace(name = :a)
  ```
- **Use keyword shorthand when forwarding same-named locals.** Write
  `f(args...; kwarg1, kwarg2)`, not
  `f(args...; kwarg1 = kwarg1, kwarg2 = kwarg2)`.

### Type system

- **No abstract-typed fields.** Every struct field must be concretely typed.

### Performance

- **Type stability first.** All operations should be inferable — verify with `@inferred`.
- **Minimize allocations in the core paths.** Pay particular attention to Fourier lowering, recursive expansion, and Liouvillian term composition. Sink type-parameters (`F` in `where {F}`) force specialization so nested `do`-blocks inline with zero closure allocation.
- **No kwargs in hot paths.** Keyword arguments prevent specialization and can allocate. Expose kwargs at the API boundary, forward to positional-arg inner functions.

### Imports and style

- **No `using X` without explicit imports.** Use `using X: func1, func2` or `import X`. ExplicitImports.jl enforces this.
- **Format with JuliaFormatter.** Run `make format` and respect
  `.JuliaFormatter.toml` before committing.

### Documentation and comments

- **Comments: compact, why not what.** Default to no comment; add one only for a non-obvious *why*. Keep it to a couple of lines.
- **No meta comments.** Never comment to explain what the code does — the code itself is the explanation. Only comment for a non-obvious *why*.
- **No docstrings on internal helpers.** Let the function name and signature speak. Docstrings belong on the public interface.

### Naming

- **No underscore-prefixed names.** Do not use `_filename.jl`, `_function_name`, or any leading-underscore convention for private internals. Julia's module system handles visibility.
- **Prefer untyped parameters for higher-order functions.** Julia's pass-through heuristic skips specialization for `Function`-annotated args that aren't called directly — but the same heuristic applies to untyped `f` too. Use `f::F where F` only when you need to force specialization (e.g. nested closures). In practice the difference is negligible; keep signatures simple.
  
## Performance terminology

Three distinct axes; keep them separate when reporting numbers.

- **Precompilation time**: compiling package code into the on-disk cache. Paid once per version/dependency change, not per session; where `precompile.jl`'s `@compile_workload` runs.
- **TTFX**: cold-process latency to the *first* call: `using` load time plus any JIT not cached during precompilation. A bigger `@compile_workload` trades more precompilation time for less TTFX.
- **Runtime**: steady-state cost of every call after the first, no compilation. What `@btime`/`benchmark/` measure after warmup; target of the zero-allocation passes.

## Agent skills

### Issue tracker

Issues and specs live in GitHub Issues for `oameye/FloquetExpansions.jl`; use the `gh` CLI. See `docs/agents/issue-tracker.md`.

### Triage labels

Use `needs-triage`, `needs-info`, `ready-to-implement`, `needs-maintainer`, and `wontfix`. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context layout with root `CONTEXT.md` and `docs/adr/`. See `docs/agents/domain.md`.

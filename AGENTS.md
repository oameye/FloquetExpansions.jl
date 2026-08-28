# FloquetExpansions.jl

## What is this?

A Julia package for the computation of FloquetExpansions included the Van-Vleck high frequency expansion and the Floquet-Magnus expansion.
 
## Git policy

**Never commit or push.** Neither Claude nor any subagent may run `git commit`, `git push`, or any git command that modifies history. All commits are made by the user, except explicity told otherwise.

## Development workflow

Common tasks via the Makefile:

```sh
make test       # run all tests (Pkg.test)
make format     # format with Runic (src/ test/ benchmark/)
make docs       # build documentation
make servedocs  # serve docs locally with LiveServer
make bench      # run benchmarks
make all        # setup + format + test + docs
```

### Quick debugging with Julia MCP

Use the `julia-mcp` MCP server (tools: `julia_eval`, `julia_list_sessions`, `julia_restart`) for quick debugging and testing small snippets — e.g., checking a type, evaluating an expression, or verifying a method signature. Prefer this over spinning up a full test run when you just need a quick answer.

### Test structure

Tests are organized in subdirectories matching `src/`, run via `test/runtests.jl` (ParallelTestRunner).

### Testing patterns

- `@inferred` for type stability checks
- `@allocations` for zero-allocation verification on hot paths

### Quality gates

Before merging any PR:
1. `make test` passes (all quality + unit + integration tests)
2. JET reports zero issues
3. No `Any`-typed fields in any struct (CheckConcreteStructs)
4. All imports explicit (ExplicitImports)

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

### Type system

- **No abstract-typed fields.** Every struct field must be concretely typed.
- **`ProductSpace{T}` uses concrete `Tuple` storage.** The type parameter `T` is a concrete tuple type.

### Performance

- **Type stability first.** All operations should be inferable — verify with `@inferred`.
- **Minimize allocations in the pipeline.** The streaming passes in `passes.jl` are the hot path. Sink type-parameters (`F` in `where {F}`) force specialization so nested `do`-blocks inline with zero closure allocation.
- **No kwargs in hot paths.** Keyword arguments prevent specialization and can allocate. Expose kwargs at the API boundary, forward to positional-arg inner functions.

### Imports and style

- **No `using X` without explicit imports.** Use `using X: func1, func2` or `import X`. ExplicitImports.jl enforces this.
- **Format with Runic.** Run `make format` before committing.

### Documentation and comments

- **Comments: compact, why not what.** Default to no comment; add one only for a non-obvious *why*. Keep it to a couple of lines.
  
## Performance terminology

Three distinct axes; keep them separate when reporting numbers.

- **Precompilation time**: compiling package code into the on-disk cache. Paid once per version/dependency change, not per session; where `precompile.jl`'s `@compile_workload` runs.
- **TTFX**: cold-process latency to the *first* call: `using` load time plus any JIT not cached during precompilation. A bigger `@compile_workload` trades more precompilation time for less TTFX.
- **Runtime**: steady-state cost of every call after the first, no compilation. What `@btime`/`benchmark/` measure after warmup; target of the zero-allocation passes.


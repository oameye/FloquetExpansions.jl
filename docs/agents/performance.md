# Performance

Authority for performance rules and for the vocabulary used when reporting numbers. `STANDARDS.md` routes here.

## The three axes

Keep these separate. A number reported without saying which axis it measures is unreadable.

- **Precompilation time**: compiling package code into the on-disk cache. Paid once per version or dependency change, not per session. Where `precompile.jl`'s `@compile_workload` runs.
- **TTFX**: cold-process latency to the *first* call, meaning `using` load time plus any JIT that precompilation did not cache. A bigger `@compile_workload` trades more precompilation time for less TTFX.
- **Runtime**: steady-state cost of every call after the first, with no compilation. What `@btime` and `benchmark/` measure after warmup, and the target of the zero-allocation passes.

Call the `julia-tffx` skill for the first two, `optimize-julia-code` and `profile-performance` for the third.

## Rules

- **Keep every operation inferable.** Assert it with `@inferred` in the test for the behaviour. JET covers the package-wide case.
- **Keep the core paths allocation-free**, and assert it with `@allocations`. The paths that matter are Fourier lowering, the recursive expansion, and Liouvillian term composition. Sinking a type parameter (`F` in `where {F}`) forces specialization, so a nested `do` block inlines with no closure allocation.
- **Take keyword arguments at the API boundary and forward to a positional inner function.** Keyword arguments prevent specialization and can allocate, so a hot path takes positional arguments only.

## Measuring

`make bench` runs the suite. `.github/workflows/Benchmarks.yaml` tracks it and comments on a pull request, alerting above 130% and failing above 170% of the recorded baseline.

Measure rather than reason about allocations or inference. A claim about either needs the number that produced it.

# Performance

Authority for inference, optimizer-stability, lowering, allocation, and runtime-performance rules, and for the vocabulary used when reporting numbers. `STANDARDS.md` routes here.

## The three axes

Keep these separate. A number reported without saying which axis it measures is ambiguous.

- **Precompilation time**: time spent compiling package code into reusable cache artifacts after a package or dependency change. The repository currently has no dedicated precompile workload or precompilation gate.
- **TTFX**: cold-process latency to the first relevant call: package loading plus JIT work not already cached. The repository currently has no TTFX gate.
- **Runtime**: steady-state execution cost after compilation. `benchmark/` and `Benchmarks.yaml` are the repository's current automated runtime-performance surface.

Do not use a runtime benchmark as evidence about TTFX, or a successful precompile as evidence about steady-state runtime.

## Inference and optimizer stability

- **Treat inference as an acceptance property on core paths.** Use `@inferred` for stable observable return types and package-wide JET for broader inference diagnostics.
- **Keep the explicit-frame completion path compiler-clean.** `test/quality/JET.jl` contains required `JET.@test_opt` workloads for explicit-frame Gram completion, recursive Gram completion, Spectral completion, and completed-state accessors. Do not remove or weaken those workloads to hide an inference problem.
- **Keep the explicit-frame computational core concrete.** Automatic dissipative-frame discovery is allowed to be a dynamic convenience frontend because the number of independent directions is discovered at runtime; once a `DissipativeFrame` is explicit, the downstream completion path should preserve concrete dispatch and storage.
- **Do not introduce boxed closure captures in `src/`.** CodeRatchet's `boxes` metric holds the current baseline; a `Core.Box` is a lowering-level sign that captured state escaped the intended concrete local representation.
- **Treat undismissed JETLS lowering diagnostics as regressions.** CodeRatchet's `lsp` metric ratchets them for `src/`. A dismissal belongs in `code_ratchet/rulings.toml` with a concrete reason; it is not a substitute for fixing a real lowering problem.

## Runtime and allocations

- **Minimize avoidable allocations; do not claim whole symbolic algorithms are allocation-free.** Fourier lowering, recursive expansion, Liouvillian composition, symbolic matrix/series algebra, and positive completion naturally construct symbolic objects. Optimize repeated temporary structure and data movement, and measure the path that matters.
- **Reuse structural work.** When repeated symbolic solves share the same leading matrix, build and reuse the solve plan. In Gram recursion, dark-sector dressing and the associated Feshbach residual should reuse equivalent solve structure rather than factor the same system twice.
- **Preserve Hermitian structure in completion linear algebra.** Use Hermitian/congruence elimination directly and avoid materializing dense elementary transforms when structured elimination suffices.
- **Keep API convenience separate from hot kernels.** Keyword arguments are valid at public boundaries. Do not mechanically ban them; use positional inner kernels where measurements or inference show that doing so keeps a compiler-sensitive call chain concrete and simple.
- **Do not change the completion scalar backend as an incidental optimization.** The current dedicated completion scalar representation remains behind `completion_conversion.jl`; replacing it with a native `SQA.CNum`-based layer is a separate architectural change.

## Measuring

`make bench` runs the runtime benchmark suite. `.github/workflows/Benchmarks.yaml` records benchmark history, comments when a result exceeds 130% of baseline, and fails when a result exceeds 170%.

`make ratchet` additionally checks source-level non-regression metrics, including boxed captures, JETLS diagnostics, and its local JET metric. Ratchet CI omits the duplicate JET metric because `JET.yml` already requires the package's absolute-zero JET gate.

A performance claim must identify the axis, workload, and measurement that supports it. Allocation claims need an allocation measurement; the repository currently has no package-wide committed `@allocations` gate, so do not imply one exists.

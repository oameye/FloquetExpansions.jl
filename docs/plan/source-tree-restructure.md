# Source tree restructure

Status: implemented in the working tree, uncommitted. Behaviour, public names, and the `@public`
seam are unchanged.

## Goal

Replace the flat `src/` (20 files, one subfolder) with layer folders that follow the data flow
in [`docs/agents/architecture.md`](../agents/architecture.md), and leave an obvious home for the
work already queued on the tracker. The test tree mirrors it.

The restructure is not a refactor of the algorithms. A file is moved, or split along a boundary
that already existed in it, and nothing else.

## Problems with the old layout

1. **Flat directory with prefix-namespacing.** Eight files carried a `completion_` or `gram_`
   prefix to do the job of a folder. `bloch_feshbach/` was the only subfolder.
2. **Include order was not layered.** Completion files were interleaved with core files.
3. **Three backward type dependencies forced that order.**
   - `liouvillian.jl` built `MicroscopicProvenance` and its seed types, which lived in
     `completion_types.jl`.
   - `FloquetExpansion{…,C<:Completion,R<:FloquetProvenance}` needed completion types first, and
     `floquet_expansion.jl` defined `effective_generator` for `C<:PositiveCompletion`.
   - `periodic_operator.jl` defined `antiderivative(G, ::VanVleck)`, so gauges had to precede
     the Fourier layer.
4. **Common completed-state accessors lived in the Gram file** (`gram_completion.jl`), while
   the architecture map assigned them to `completion.jl` and Spectral relies on them too.
5. **Mixed files.** `liouvillian.jl` held channels, displays, rate validation, provenance, and
   the map algebra. `gksl_coordinates.jl` held the frame type with its linear algebra and the
   GKSL extraction.
6. **Names that no longer matched.** `periodic_operator.jl` defined `PeriodicGenerator`;
   `test/engine.jl` tested Hori–Deprit and `test/collector.jl` tested Fourier lowering.

## Layout

One module, plain folders, no Julia submodules. Submodules would add qualified-name plumbing and
ExplicitImports churn for no user-visible gain.

```text
src/
├── FloquetExpansions.jl
├── generators/                  periodic generators; no gauges, expansions, or completion
│   ├── periodic_generator.jl    ← periodic_operator.jl, minus the gauge-fixed antiderivative
│   ├── liouvillian.jl           ← liouvillian.jl: map algebra, compose, Fourier lowering
│   ├── channels.jl              ← liouvillian.jl: channels, displays, rate validation,
│   │                               provenance construction; + provenance types
│   └── quasienergy.jl
├── words/                       harmonic word algebra, generic in its coefficient type
│   ├── lyndon_words.jl          ← bloch_feshbach/lyndon_words.jl
│   └── van_vleck_words.jl       ← bloch_feshbach/van_vleck_words.jl
├── expansion/
│   ├── gauges.jl                + antiderivative(G, ::VanVleck)
│   ├── floquet_expansion.jl     + Completion, Uncompleted
│   ├── hori_deprit.jl
│   └── bloch_feshbach.jl        ← bloch_feshbach/bloch_feshbach.jl
├── gksl/
│   ├── dissipative_frame.jl     ← gksl_coordinates.jl lines 1-335
│   ├── coordinates.jl           ← gksl_coordinates.jl lines 337-605
│   └── floquet.jl               ← gksl_floquet.jl, minus the completed hamiltonian method
└── completion/
    ├── types.jl                 ← completion_types.jl, minus state supertype and provenance
    ├── backend/                 the dedicated completion scalar backend (#85 replaces this)
    │   ├── matrix_series.jl
    │   ├── linear_algebra.jl    ← completion_linear_algebra.jl
    │   └── conversion.jl        ← completion_conversion.jl
    ├── frame_discovery.jl       ← completion_frame.jl
    ├── accessors.jl             new: common completed accessors, stored_* readers,
    │                               completed effective_generator and hamiltonian
    ├── gram/
    │   ├── factorization.jl     ← gram_completion.jl, minus the common accessors
    │   └── recursion.jl         ← gram_recursion.jl
    ├── spectral.jl              ← spectral_completion.jl
    └── positive_completion.jl   ← completion.jl, minus the stored_* readers
```

Tests mirror it: `test/generators/{periodic_generator,fourier_lowering,liouvillian,quasienergy}`,
`test/expansion/{hori_deprit,expansion_algorithms,liouvillian_phase,kpo}`, `test/gksl/coordinates`,
and `test/completion/{state,storage,gram,gram_recursion,spectral,residual_scaling,validation,
analytic_examples,numeric_oracle}`.

The layering rule and the per-file ownership table live in `docs/agents/architecture.md`, which
is their authority.

## Why this shape: the queued work

- **Shared projection core** (#102, #182, #115, and the Bloch projection files in draft #324).
  The word algebra is already generic over the coefficient type and never touches a generator
  component, so it becomes `words/`. The next consumer, CK output amplitudes or QHB, imports
  it rather than reaching into the Bloch/Feshbach algorithm.
- **Native GKSL expansion** (#338 and its stack #337, #339 to #342). It builds a GKSL generator
  directly and deliberately never calls `positive_completion`. It gets its own folder after
  `gksl/`, with no dependency on `completion/`.
- **CK period map** (draft #324). Map-level reconstruction is a separate operation from the
  static expansion. It gets its own folder when its public API exists.
- **Floquet Schrieffer–Wolff** (#311) and further algorithms. One file in `expansion/` plus a
  selector or gauge in `expansion/gauges.jl`. Moving the gauge-fixed `antiderivative` there
  means a new gauge adds its integration constant beside its type.
- **Native CNum completion scalars** (#85). The scalar backend is isolated in
  `completion/backend/`, so that change replaces one folder.
- **CP-completed default** (#313). `floquet_expansion` would call `positive_completion`. That
  is a run-time call into a later folder, which the layering rule allows. No type moves.
- **Multi-frequency harmonics** (#107, #108). They change `generators/periodic_generator.jl` and
  the letters in `words/`. Everything above sees harmonics only through `PeriodicGenerator`.

## Verification done

- Before any edit, all 840 names and 566 method signatures defined by the module were recorded.
  After the restructure every named binding and every signature is identical. Only
  compiler-generated closure names differ, because their numbering follows definition order, and
  their count is unchanged.
- `make test`: 898 of 898 pass, with all 23 discovered test files running under their new paths.
- `make jet`: the package JET report and both optimizer-stability testsets pass.
- `make format`: no changes.
- CodeRatchet `all refresh` passes all five gates and regenerated the path-keyed baselines.
  The complexity sums (`cyc_sum` 950, `cog_sum` 645, `arg_sum` 740) equal a fresh measurement
  of `main` at 8fe4108, so the split moved code without changing it. The committed baseline had
  been stale at 954, 647 and 752. CodeRatchet refuses a git worktree because it checks for a
  `.git` directory, so the refresh ran in a scratch clone of this tree.

## Still to do by the maintainer

- `make docs`. No `@autodocs Pages=` filter exists, so the manual should be unaffected, but only
  the build proves it.
- Commit. Suggested split: the `git mv` renames, then the relocations and splits, then the docs
  and ratchet baselines.

## Rebasing the open PRs

Every file moved, so the native GKSL stack and #336 conflict. `git rebase` with rename detection
carries edits to moved files over. New files need placing by hand:

| PR | New file | Goes to |
| --- | --- | --- |
| #338 | `src/native_gksl.jl` | `src/native_gksl/` folder after `gksl/` in the include list |
| #338 | `test/native_gksl_*.jl` | `test/native_gksl/` |
| #337, #339 to #342 | `test/native_gksl_*.jl` | `test/native_gksl/` |
| #336 | edits `src/completion_types.jl` | `src/completion/types.jl` |
| #342 | edits `test/cp_completion_validation.jl` | `test/completion/validation.jl` |

After each rebase, refresh the ratchet baselines again.

## Decisions taken

1. Internal identifiers keep their names, including now-redundant prefixes such as
   `completion_matrix_zeros`. The prefix still disambiguates inside a single module.
2. ADR verification pointers to test files were amended to the new paths, with no other change
   to any ADR.
3. `completion/backend/matrix_series.jl` (613 lines) stays whole. It is coherent, and #85 will
   rewrite it anyway.

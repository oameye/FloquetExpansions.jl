# van Vleck expansion — working notes

| file | what it is | read it when |
|---|---|---|
| `DESIGN.md` | The design. Types, operations, collector, engine, API, tests, performance policy. | You are implementing. |
| `DERIVATION.md` | The recursion derived from the spec, and its verification through order 5. | You need the mathematics. |
| `DECISIONS.md` | Decision log D1–D7, plus paths evaluated and rejected, with reasons. | Before re-proposing anything. |
| `RESEARCH.md` | External findings: SQA API surface, the published oracle, adversarial review, algorithms and tooling, gauge experiments. | You need a fact and its source. |

Runnable experiments (isolated env, not part of the package):

| script | what it establishes |
|---|---|
| `collector_prototype.jl` | Harmonic extraction from a symbolic `QAdd`, round-trip exact. |
| `gauge_experiment.jl` | van Vleck vs des Cloizeaux; they part at order 3. |
| `nikolaev_experiment.jl` | Nikolaev's explicit formula: orders 0–2 exact, order 3 unresolved. |

Source papers are NOT tracked; they live in `tmp/` locally and are all public:

| file in `tmp/` | what | get it from |
|---|---|---|
| `floquet_bath_tracing.tex` | the authoritative spec, §4 (`sec:vv`) | the author |
| `2108.02861v2.pdf` | oracle tables; conventions must be translated | arXiv:2108.02861 |
| `nikolaev_1612.05207.pdf` | parked fast-path candidate | arXiv:1612.05207 |

One rule that cost several wrong turns to learn, and that governs any future engine work:

> `eq:Heff1` and `eq:Heff2` CANNOT discriminate gauges. Every canonical scheme tested agrees with
> van Vleck through order 2 and parts company at order 3. Validate at order 3 or beyond, and gate
> on ISOSPECTRALITY first: `Heff`'s spectrum is the physical quasienergy spectrum and is gauge
> invariant, so spectra differing means the implementation is wrong, while spectra matching and
> operators differing means a gauge difference (which converts cheaply).

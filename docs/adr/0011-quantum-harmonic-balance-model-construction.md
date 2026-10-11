# Quantum harmonic balance is a model construction in its own layer

Quantum harmonic balance (QHB) lives in `src/harmonic_balance/`, after `gksl/`. It maps a physical Hamiltonian or `Liouvillian` of one bosonic mode to a time-independent generator on carrier modes, plus a reconstruction map back to physical observables. It is not a gauge, not an expansion algorithm, and not a `FloquetExpansion`.

- **Not a gauge or an algorithm.** ADR 0010 requires every algorithm for a gauge to return the same retained coefficients, and a gauge fixes one Floquet effective generator. QHB acts on an enlarged algebra, so neither applies.
- **Not a `FloquetExpansion`.** ADR 0004 keeps one result model for high-frequency expansions. QHB returns an ordinary Hamiltonian or `Liouvillian`, and its reconstruction lives on the `CarrierEmbedding`, so ADR 0004 is unchanged. Reconstruction is not micromotion.
- **Its own folder.** The architecture map sends an open-system construction that produces a generator directly to a folder after `gksl/` that does not depend on `completion/`.

## The construction

Every carrier frequency ``ω_j`` contributes a co-rotating and a counter-rotating ladder mode:

```math
I_t(a) = λ_H^{-1/2} \sum_j \left(e^{-iω_j t} a_{j+} + e^{+iω_j t} a_{j-}\right),
\qquad λ_H = 2N_h.
```

- **Why it is exact.** The embedding is a unitary change of frame on the physical mode, extended by dark partner modes. So it preserves every product of the physical algebra, provided the weight ``λ_H^{-1/2}`` is exact. Exact weights need SecondQuantizedAlgebra's exact numeric radicals (qojulia/SecondQuantizedAlgebra.jl#284).
- **The generator.** The carrier generator is ``λ_H P[I_t(\cdot)] + H_{\mathrm{fr}}``:
  - ``P`` keeps the processes whose frequency vanishes identically;
  - ``H_{\mathrm{fr}} = -\sum_j ω_j (n_{j+} - n_{j-})`` is the gauge term of the rotating frame;
  - for a `Liouvillian`, ``P`` acts on each ``ρ ↦ AρB`` term as a whole.
- **Why ``λ_H`` is required.** The dark partners never reach physical observables, so their Hamiltonian is a gauge choice. With a zero dark Hamiltonian, the rotating-frame generator keeps terms of the size of the physical frequency, and averaging them away gives the wrong spectrum. Giving every dark partner a copy of the physical generator makes the rotating-frame generator exactly time independent for quadratic Hamiltonians and linear dissipators, and that copy is what ``λ_H P[I_t(\cdot)]`` produces. For nonlinear terms the copies do not commute with the carrier rotation, and ``P`` is an approximation.
- **Frequencies.** Carrier and drive frequencies are symbolic expressions with exact coefficients, and relations between carriers are written into those expressions, such as ``(ω_1 + ω_2)/2``. A resonance is a frequency that expands to zero. Floating-point coefficients of a symbol are rejected, because they make that zero test unreliable.
- **Scope.** The first implementation embeds one bosonic mode and rejects generators with any other operator. The reduction to QHB-RWA and higher-order corrections are later work.

## Considered options

- **A quadrature carrier basis** ``(u_j, v_j, p_{u_j}, p_{v_j})`` with the rescaled commutator ``[u, p_u] = iħ/λ_H``, as in the derivation notes. Rejected: the ladder modes ``a_{j±}`` are a passive change of basis of the same variables, and they keep the standard commutators of SecondQuantizedAlgebra.
- **One ladder mode per carrier.** Rejected as the primary construction, because it keeps only the co-rotating sector. It is the later QHB-RWA reduction, not unreduced QHB.
- **Embedding plus a plain time average, without ``λ_H``.** Rejected, because it gives the wrong spectrum (see above).
- **A `resonances` keyword** listing retained frequency combinations. Rejected, because relations written into the carrier expressions determine the resonances without a second declaration.

## Gate

`make test` runs `test/harmonic_balance/unreduced.jl`, whose testsets hold this decision:

- "the carrier embedding preserves the physical algebra"
- "number-conserving linear dynamics reconstructs exactly"
- "damped counter-rotating dynamics keeps the physical poles"
- "resonances follow the declared carrier frequencies"
- "invalid embeddings and generators are rejected"

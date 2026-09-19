# Native CP-HFE internal development note

This page is intentionally not part of the public manual navigation. It records the internal data flow exercised by issue #105 while the public API remains undecided.

The first implementation tranche keeps physical dissipative amplitudes separate from generic Liouvillian algebra:

```text
physical collapse/jump channel
    -> periodic physical amplitude
    -> coherent Van Vleck kick transport
    -> retained amplitude series
    -> finite transported amplitude
    -> one Kraus row per Fourier harmonic
    -> direct GKSL generator
```

For a physical amplitude `L_a(t)` transported by the coherent kick,

```math
L'_a(t)=e^{iK(t)}L_a(t)e^{-iK(t)},
```

the finite retained amplitude is squared only after truncation. Period averaging then gives

```math
\overline{\mathcal D[L'_a(t)]}
=\sum_m \mathcal D[L'_{a,m}],
```

so the finite generator is manifestly GKSL without constructing or repairing a truncated generic-Liouvillian Kossakowski matrix.

This is distinct from `positive_completion`: the amplitude series is primary perturbative data and the final finite Gram/Kraus square is not truncated again. The current tranche accepts periodic collapse/jump operators with time-independent external jump rates; periodic scalar rates remain outside scope until their amplitude-level square-root/onset semantics are specified.

No user-facing selector or result type is committed by this tranche. The implementation remains an internal consumer used to settle the #72 architecture seam.

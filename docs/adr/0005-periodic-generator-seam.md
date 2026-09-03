# `PeriodicGenerator` is the Fourier boundary

`PeriodicGenerator` is the canonical representation of a periodic generator: it carries the drive frequency and its nonzero integer harmonics, while missing harmonics evaluate as zero. It replaces the Hamiltonian-only `PeriodicOperator`; no alias or separate periodic Liouvillian container is maintained.

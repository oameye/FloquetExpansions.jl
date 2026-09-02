# Quasienergy operators include dissipative generators

`QuasienergyOperator` accepts both Hamiltonian and Liouvillian periodic generators. For a
Hamiltonian it stores the Sambe operator with blocks `H[m - n] - m*wd*δ[m, n]`. For a Liouvillian
it stores the energy-like Floquet-Liouville operator with blocks
`im*L[m - n] - m*wd*δ[m, n]`. The latter produces generally complex dissipative quasienergies;
their imaginary parts describe decay or growth.

The type stores symbolic blocks only. Numerical vectorization and eigenspectrum computation are
separate adapters and are not part of this module.

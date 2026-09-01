# Raw arbitrary-order effective generator

The package returns the direct arbitrary-order Ikeda van Vleck expansion as the effective generator and its micromotion. This is the faithful algebraic result for Hamiltonian and Lindblad inputs; a finite-order effective Liouvillian is not promised to retain GKSL form or complete positivity.

Positivity repair, Kossakowski reconstruction, symbolic inequalities, and certificates are outside the initial implementation scope. The package must not silently replace the raw expansion with a different completely-positive approximation.

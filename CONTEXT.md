# FloquetExpansions.jl

Shared vocabulary for symbolic high-frequency expansions of periodically driven quantum systems.

## Dynamics

**Periodically driven Lindblad system**:
An open quantum system whose Lindblad generator varies periodically with time.

**Lindblad generator**:
A density-operator generator written as a Hamiltonian action plus dissipative channels with physically valid rates.
_Avoid_: using this term for a finite-order effective result unless its Lindblad form is known.

**Liouvillian**:
A linear map that generates or approximates density-operator evolution; it need not have Lindblad form.
_Avoid_: treating every Liouvillian as completely positive.

**Hamiltonian action**:
The coherent contribution ``-i[H, ρ]`` to a density-operator generator.

**Dissipator**:
The map ``D[L](ρ) = LρL† - (L†Lρ + ρL†L)/2`` associated with a collapse operator.

**Collapse operator**:
An operator used directly in ``D[L]``; any amplitude or phase belonging to that channel is included in ``L``.

**Rate-weighted jump channel**:
A bare jump operator and a separate scalar rate contributing ``γD[J]``; the rate is not part of ``J``.

**Time-dependent dissipative channel**:
A dissipative channel whose operator or scalar rate varies periodically with the drive.

## Floquet expansion

**Periodic generator**:
A time-periodic generator expressed through Fourier harmonics in a drive-frequency basis.
_Avoid_: `PeriodicOperator`, the superseded Hamiltonian-only name.

**Van Vleck expansion**:
A high-frequency expansion that separates a periodic generator into a time-independent effective generator and periodic micromotion in a chosen gauge.

**Floquet expansion**:
The finite-order result of applying a high-frequency expansion to a periodic generator, including its effective generator and micromotion.

**Effective generator**:
The time-independent generator that approximates the slow dynamics after the expansion.

**Finite-order effective generator**:
A truncated algebraic result that may not retain Lindblad form or complete positivity.

**Micromotion**:
The periodic transformation that relates the original time-dependent dynamics to the effective dynamics; its generator need not be unitary for dissipative systems.
_Avoid_: calling dissipative micromotion a unitary kick operator.

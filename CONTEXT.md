# FloquetExpansions.jl

Domain vocabulary for symbolic high-frequency expansions of periodically driven open quantum systems.

## Dynamics

**Periodically driven Lindblad system**:
An open quantum system whose density operator evolves under a Markovian Lindblad generator that is periodic in time.
_Avoid_: periodic Liouvillian when the Lindblad structure is part of the model.

**Lindblad generator**:
A generator of density-operator evolution written as a Hamiltonian contribution plus dissipative channels.
The exact generator has the usual trace-preservation and complete-positivity guarantees.

**Liouvillian**:
A linear generator acting on density operators. A Liouvillian may be more general than a generator presented in Lindblad form, particularly after finite-order algebraic expansions.
_Avoid_: using Liouvillian and Lindblad generator interchangeably when a physical guarantee matters.

**Coherent generator action**:
The Hamiltonian contribution ``-i[H, \\cdot]`` to a density-operator generator.
It is a contribution to the common generator algebra, not a separate Floquet API.

**Dissipator**:
The channel action ``D[L](\\rho) = L\\rho L^\\dagger - \\frac{1}{2}\\{L^\\dagger L,\\rho\\}`` associated with a collapse operator.

**Complete collapse operator**:
An operator appearing directly in a dissipator ``D[L]``. Any amplitude folded into `L` belongs to the collapse operator.

**Rate-weighted jump channel**:
A jump operator together with a scalar rate, contributing ``\\gamma D[J]`` to the generator.
The initial symbolic API accepts ordinary scalar rates without representing or certifying sign conditions.

**Time-dependent dissipative channel**:
A dissipative channel whose operator, scalar rate, or both vary periodically in time.
Both complete collapse operators and rate-weighted jump channels are supported input forms.

## Floquet expansion

**Lindblad van Vleck expansion**:
The arbitrary-order high-frequency expansion of a periodically driven Lindblad generator, following Ikeda, Chinzei, and Sato (2021).
_Avoid_: Floquet–Markov expansion, which implies a particular microscopic bath construction.

**Periodic generator**:
A time-periodic dynamical generator represented by its Fourier components in the drive-frequency basis.
The same container is used for Hamiltonian and Liouvillian components.
_Avoid_: `PeriodicOperator`, which is the superseded Hamiltonian-only name.

**Floquet expansion**:
An order-resolved van Vleck decomposition of a periodic generator into effective generator components and periodic micromotion.

**Effective generator**:
The time-independent generator produced by the van Vleck expansion.
For a Hamiltonian input it is an effective Hamiltonian action; for a Lindblad input it is an effective Liouvillian.

**Raw effective generator**:
The direct arbitrary-order van Vleck result. A finite-order result is not promised to retain GKSL form or complete positivity.
The initial package returns this raw result and does not apply a positivity repair.

**Micromotion**:
The periodic generator in the van Vleck decomposition that relates the original evolution to the effective generator.
_Avoid_: unitary kick operator when the generator is dissipative.

## Scope

**Symbolic-only support**:
The supported inputs and outputs are symbolic operator expressions and symbolic scalar coefficients.
Numerical matrix representations and positivity certificates are outside the initial scope.

~~~@meta
CurrentModule = FloquetExpansions
CollapsedDocStrings = true
~~~

# Floquet theory

Floquet theory is the temporal analogue of Bloch theory. It separates a periodic evolution into a
time-independent part that accumulates from cycle to cycle and a periodic part describing the
motion within one cycle. This construction is exact; the high-frequency expansion in
[High-frequency expansion](high_frequency_expansion.md) is an approximation applied afterwards.

Unless stated otherwise, the equations below use frequency units, ``\hbar=1``, as in the package's
quasienergy blocks.

## Floquet theorem

Consider a finite-dimensional, or finitely truncated, linear periodic system

~~~math
\dot{x}(t)=\mathcal{G}(t)x(t),
\qquad
\mathcal{G}(t+T)=\mathcal{G}(t),
\qquad
\omega=\frac{2\pi}{T}.
~~~

Let ``\mathcal{V}(t_2,t_1)`` be its propagator. Floquet's theorem states that there is a constant
generator ``\mathcal{G}_{\mathrm{F}}`` and a periodic, invertible map ``\mathcal{M}(t)`` such that

~~~math
\mathcal{V}(t_2,t_1)
=\mathcal{M}(t_2)
e^{(t_2-t_1)\mathcal{G}_{\mathrm{F}}}
\mathcal{M}(t_1)^{-1},
\qquad
\mathcal{M}(t+T)=\mathcal{M}(t).
~~~

Equivalently, ``\mathcal{G}_{\mathrm{F}}`` is a choice of logarithm of the one-period propagator
up to a similarity transformation. The periodic factor contains micromotion, while the constant
factor describes the stroboscopic evolution. For an open system, ``\mathcal{G}_{\mathrm{F}}`` may
be a general effective Liouvillian rather than a generator in GKLS form.

The package uses the Fourier convention

~~~math
\mathcal{G}(t)=\sum_{m\in\mathbb{Z}}\mathcal{G}_m e^{-im\omega t}.
~~~

[`PeriodicGenerator`](@ref) stores these harmonics together with the drive frequency. Its zeroth
harmonic is the period average, and missing harmonics are zero. Use [`harmonics`](@ref) to convert a
symbolic Hamiltonian or Liouvillian.

## Hamiltonian Floquet states

For a closed system, the Schrödinger equation has solutions of the form

~~~math
|\Psi_n(t)\rangle
=e^{-i\varepsilon_n t}|u_n(t)\rangle,
\qquad
|u_n(t+T)\rangle=|u_n(t)\rangle.
~~~

The periodic state ``|u_n(t)\rangle`` is the Floquet mode and ``\varepsilon_n`` is its
quasienergy. The pair is not unique: for any integer ``m``,

~~~math
\varepsilon_n\mapsto\varepsilon_n+m\omega,
\qquad
|u_n(t)\rangle\mapsto e^{im\omega t}|u_n(t)\rangle
~~~

describes the same physical state. Quasienergies are therefore defined modulo ``\omega`` and can
be restricted to a Floquet--Brillouin zone. This branch freedom is the temporal counterpart of
reciprocal-lattice equivalence in Bloch theory [Zeldovich1967, Ritus1967](@cite).

Writing the propagator in Hamiltonian form gives

~~~math
U(t_2,t_1)
=e^{-iK(t_2)}
e^{-i(t_2-t_1)H_{\mathrm{eff}}}
e^{iK(t_1)},
~~~

where ``K(t)`` is a Hermitian kick operator. The generic periodic map is instead written as
``\mathcal{M}(t)=e^{\mathcal{K}(t)}``; for dissipative dynamics its generator ``\mathcal{K}(t)``
need not be unitary or Hermitian.

## Floquet and Sambe space

Substitution of the Floquet ansatz into the Schrödinger equation gives the time-independent
quasienergy problem

~~~math
Q|u_n\rangle=\varepsilon_n|u_n\rangle,
\qquad
Q=H(t)-i\partial_t.
~~~

The operator ``Q`` acts on the tensor product of the physical Hilbert space and the space of
periodic functions, called Floquet or Sambe space. Expanding in Fourier sectors turns the periodic
problem into an infinite block matrix. With the package convention, its Hamiltonian blocks are

~~~math
Q_{mn}=H_{m-n}-m\omega_d\,\delta_{mn}.
~~~

The diagonal blocks are replicas of the averaged Hamiltonian, shifted by integer multiples of the
drive frequency; nonzero harmonics couple different replicas. [`QuasienergyOperator`](@ref)
constructs a finite symbolic truncation of these blocks. This exact time-domain-to-Sambe-space
mapping is the formulation introduced by Shirley and placed on an extended-space footing by Sambe
[Shirley1965, Sambe1973](@cite).

For a periodically driven Liouvillian ``\mathcal{L}(t)``, the corresponding energy-like operator
has blocks

~~~math
Q_{mn}=i\mathcal{L}_{m-n}-m\omega_d\,\delta_{mn}.
~~~

Its eigenvalues are generally complex: their real parts describe oscillation and their imaginary
parts describe decay or growth. They are still defined modulo ``\omega_d``. Thus ordinary
quasienergies and dissipative quasienergies are two spectral versions of the same Floquet
construction, with unitarity forcing the former to be real. The Floquet-space perspective is also
the natural starting point for the van Vleck block diagonalization discussed in
[High-frequency expansion](high_frequency_expansion.md) [Eckardt2015](@cite).

The thesis and much of the Floquet literature use ``e^{+im\omega t}`` instead. Comparing those
formulas with the package requires ``m\mapsto-m``.

## Open-system Floquet generators and complete positivity

For a periodic Markovian master equation, the propagator over each finite time interval is a
completely positive trace-preserving map. This does **not** imply that an arbitrary logarithm of the
one-period channel is a time-independent GKLS generator. The exact Floquet logarithm depends on a
logarithm branch and on the periodic similarity gauge, and the resulting effective Liouvillian need
not lie inside the GKLS cone [Schnell2021](@cite).

The same distinction appears perturbatively. A high-frequency expansion constructs a formal
Floquet Liouvillian

~~~math
\mathcal L_{\mathrm{eff}}^{[N]}
=\sum_{n=0}^{N}\omega_d^{-n}\mathcal L^{(n)},
~~~

whose retained coefficients are fixed by the microscopic periodic dynamics. Truncating this formal
series can produce a Kossakowski form with weak negative directions even though the original
time-dependent generator is Lindbladian at every time.

[`positive_completion`](@ref) addresses this finite-order problem by supplying higher-order
dissipative data that restore a positive Kossakowski form while preserving every retained
coefficient:

~~~math
\mathcal L_{\mathrm{CP}}^{[N]}
-\mathcal L_{\mathrm{eff}}^{[N]}
=\mathcal O(\omega_d^{-N-1}).
~~~

This is a perturbative completion, not an assertion that
``\mathcal L_{\mathrm{CP}}^{[N]}`` equals an exact Floquet GKLS logarithm. It also does not promote
the asymptotic high-frequency expansion into a convergent series. What is controlled is the
retained perturbative data: the completed and raw models agree through the claimed order.

Because the correction starts beyond that order, the retained periodic micromotion does not need to
be recomputed. The finite completed approximation can use the same retained kick/micromotion and a
completed effective generator; a different kick would only enter as part of higher-order Floquet
data.

### Frames, rates, and physical invariants

A Kossakowski matrix is a representation of the dissipative Hermitian form in a chosen operator
frame. Changing the [`DissipativeFrame`](@ref) changes matrix entries and can redistribute the same
physical dissipative form among different jump representatives. Likewise, a particular set of jump
operators or branch rates is not unique: unitary channel rotations, nonunitary coordinate changes,
and different positive factorizations give different representatives.

The frame-independent content is the Hermitian form itself. Positivity and inertia are invariant
under invertible congruence transformations, whereas individual rates and jump operators are
gauge/representation dependent.

This is also why different completion algorithms can be physically equivalent through retained
order but differ as finite expressions. [`Gram`](@ref) fixes one algebraic jump-amplitude gauge;
[`Spectral`](@ref) fixes a restricted perturbative decay-rate/HCM gauge. Their higher-order
continuations need not coincide.

See [Positive completion](@ref) for the user workflow and
[High-frequency expansion](high_frequency_expansion.md) for the active/dark, Feshbach, and onset
construction.

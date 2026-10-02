~~~@meta
CurrentModule = FloquetExpansions
CollapsedDocStrings = true
~~~

# [Bloch/Feshbach projection](@id bloch-feshbach-theory)

The fast Fourier sectors of a periodic problem can be eliminated by projection instead of by a
frame transformation. One projects the Sambe problem onto the zero Fourier sector and solves for a
wave operator that maps this sector onto the invariant subspace of the driven dynamics. Feshbach's
elimination of the other sectors gives an effective generator that depends on the quasienergy
[Feshbach1958, Feshbach1962](@cite), and Bloch's wave-operator equation removes that dependence
[Bloch1958](@cite). For a periodic drive this is the Brillouin--Wigner high-frequency expansion of
Mikami et al. [Mikami2016Brillouin](@cite). Its effective generator is a different representative
of the slow dynamics, and a static normalization recovers the van Vleck pair, the effective
generator and micromotion generator of the
[high-frequency expansion](@ref high-frequency-expansion-theory). Reviews of the wave-operator
formalism are [Killingbeck2003, Jolicard2003](@cite).

## Feshbach projection in Sambe space

Use the Fourier convention and the quasienergy operator of [Floquet theory](@ref floquet-theory),

~~~math
Q_{mn}=H_{m-n}-m\omega_d\,\delta_{mn}.
~~~

Let ``P`` project onto the zero Fourier sector and ``\bar P=1-P`` onto all others. The diagonal
``-m\omega_d`` makes ``P`` an exactly degenerate model space, separated from every other sector by
at least ``\omega_d``. Eliminating ``\bar P`` gives Feshbach's effective Hamiltonian on the zero sector,

~~~math
H_{\mathrm{eff}}(\varepsilon)
=H_0+\sum_{m,n\ne0}H_{-m}\left[(\varepsilon-\bar PQ\bar P)^{-1}\right]_{mn}H_n,
~~~

whose quasienergies solve the nonlinear eigenproblem
``H_{\mathrm{eff}}(\varepsilon)v=\varepsilon v``. This Brillouin--Wigner form is exact, but it
depends on the eigenvalue being sought [Feshbach1958, Mikami2016Brillouin](@cite).

## Bloch wave operator

The Bloch wave operator ``\Omega`` maps the zero sector onto the invariant subspace connected to it,
with Bloch's intermediate normalization ``P\Omega P=P``. It satisfies

~~~math
Q\,\Omega=\Omega\,PQ\,\Omega,
~~~

which replaces the eigenvalue ``\varepsilon`` of the Feshbach problem by right multiplication with
an effective operator and is therefore energy independent [Bloch1958, Mikami2016Brillouin](@cite).
The ``(m,0)`` blocks of ``\Omega`` are the Fourier components of a periodic map ``\Omega(t)``. In
the generic map notation ``\dot x=\mathcal{G}(t)x``, with ``\mathcal{G}=-iH`` for a Hamiltonian,
the equation reads

~~~math
\partial_t\Omega=\mathcal{G}\,\Omega-\Omega\,\mathcal{B},
\qquad
\langle\Omega\rangle=1,
\qquad
\mathcal{B}=\langle\mathcal{G}\,\Omega\rangle.
~~~

Solutions of the slow problem ``\dot y=\mathcal{B}y`` map to driven solutions ``x=\Omega y``. We
call ``\mathcal{B}`` the Bloch effective generator; for a Hamiltonian, ``H_{\mathrm{B}}=i\mathcal{B}``
is the energy-independent Brillouin--Wigner Hamiltonian of Mikami et al. Substituting
``\mathcal{B}=\langle\mathcal{G}\,\Omega\rangle`` makes the equation quadratic in ``\Omega``.

## Recurrence

Expand in inverse powers of the drive frequency,

~~~math
\Omega=\sum_{n\ge0}\omega^{-n}\Omega^{(n)},
\qquad
\Omega^{(0)}=1,
\qquad
\langle\Omega^{(n)}\rangle=0\quad(n\ge1),
\qquad
\mathcal{B}=\sum_{n\ge0}\omega^{-n}\mathcal{B}^{(n)}.
~~~

The Bloch equation becomes

~~~math
\mathcal{R}^{(n)}
=\mathcal{G}\,\Omega^{(n)}
-\sum_{j=1}^{n}\Omega^{(j)}\mathcal{B}^{(n-j)},
\qquad
\mathcal{B}^{(n)}=\langle\mathcal{R}^{(n)}\rangle,
\qquad
\Omega^{(n+1)}_m=\frac{i}{m}\mathcal{R}^{(n)}_m\quad(m\ne0).
~~~

The last step is the zero-average inverse of ``\partial_\tau`` from the homological equation. Each
order is a linear step whose residual collects products of known lower-order coefficients; no
nested commutators with a transformation generator appear. Mikami et al. derived the same
recursion for Hamiltonians in their Floquet Brillouin--Wigner expansion
[Mikami2016Brillouin](@cite). In the static setting,
Bravyi, DiVincenzo, and Loss note that each coefficient of Bloch's series is a second-degree
polynomial in the lower-order ones, a simpler structure than the Schrieffer--Wolff series
[Bravyi2011](@cite).

## Leading terms and non-Hermiticity

The first two coefficients equal the van Vleck ones,

~~~math
\mathcal{B}^{(0)}=\mathcal{G}_0,
\qquad
\mathcal{B}^{(1)}=i\sum_{m>0}\frac{[\mathcal{G}_{-m},\mathcal{G}_m]}{m}.
~~~

From order ``\omega^{-2}``, ``\mathcal{B}`` is a different representative with the same
quasienergies. For a Hamiltonian, ``H_{\mathrm{B}}`` is then in general not Hermitian: the
zero-sector projections of the Floquet states are not orthonormal, and their metric is
``S=\langle\Omega^\dagger\Omega\rangle`` [Bloch1958, Mikami2016Brillouin, Kvaal2008](@cite).
des Cloizeaux's similarity ``S^{1/2}H_{\mathrm{B}}S^{-1/2}`` restores Hermiticity. For the split
of the zero sector against all others it coincides with the Schrieffer--Wolff effective
Hamiltonian and with Shavitt and Redmon's canonical form of static quasidegenerate perturbation
theory applied to that two-block split [DesCloizeaux1960, ShavittRedmon1980, Kvaal2008, Bravyi2011](@cite).

## Normalization to the van Vleck gauge

The wave operator and the van Vleck micromotion describe the same invariant subspace in different
normalizations. Seek a constant ``N``, the static factor, and a periodic ``\Lambda``, the connected
logarithm, with

~~~math
\Omega(t)N=e^{\Lambda(t)},
\qquad
\langle\Lambda\rangle=0.
~~~

Averaging with ``\langle\Omega\rangle=1`` fixes ``N=\langle e^{\Lambda}\rangle``. Substituting
``\Omega=e^{\Lambda}N^{-1}`` into the Bloch equation gives

~~~math
e^{-\Lambda}\left(\mathcal{G}\,e^{\Lambda}-\partial_t e^{\Lambda}\right)
=N^{-1}\mathcal{B}N.
~~~

The left side is the Lie transform generated by ``\Lambda``, and ``\langle\Lambda\rangle=0`` is the
van Vleck condition. Uniqueness of the van Vleck pair therefore identifies

~~~math
\mathcal{K}=\Lambda,
\qquad
\mathcal{G}_{\mathrm{eff}}=N^{-1}\mathcal{B}N.
~~~

The relations Mikami et al. give between their Brillouin--Wigner, Floquet--Magnus, and van Vleck
Hamiltonians imply the same static factor [Mikami2016Brillouin](@cite). It starts at second order,

~~~math
N=1+\frac{1}{2\omega^{2}}\sum_{m\ne0}\frac{\mathcal{G}_m\mathcal{G}_{-m}}{m^2}
+\mathcal{O}(\omega^{-3}),
\qquad
\mathcal{G}_{\mathrm{eff}}^{(2)}=\mathcal{B}^{(2)}+\left[\mathcal{G}_0,N^{(2)}\right].
~~~

Since ``\Lambda`` equals the Lie-transform micromotion, it is a Lie series in the harmonics
``\mathcal{G}_m``: every term is a nested commutator, and all other products cancel. The same
holds for ``\mathcal{G}_{\mathrm{eff}}=N^{-1}\mathcal{B}N``, which equals the Lie-transform
effective generator, although ``\mathcal{B}`` and ``N`` separately are not Lie series.

For a Hamiltonian, ``N=S^{-1/2}V`` with ``V`` unitary, so the van Vleck and des Cloizeaux
Hamiltonians differ by a static unitary rotation of the zero sector. They are not equal in general.
The van Vleck gauge decouples every pair of Fourier sectors, which is the canonical form applied
to all sectors at once [Eckardt2015](@cite), whereas des Cloizeaux's normalization only decouples
the zero sector from the rest. The first difference appears at order ``\omega^{-3}`` or later, depending on which harmonics
the drive contains.

## Open systems and complete positivity

For a Liouvillian, ``\mathcal{G}=\mathcal{L}``, and every step above is algebraic: ``\Omega``,
``N``, and ``\Lambda`` are general maps, and ``N^{-1}\mathcal{B}N`` is a similarity rather than a
unitary rotation. The split into the zero sector and the rest is then a projection-operator
technique in the sense of Nakajima and Zwanzig, with the period average as the relevant part, and
Feshbach's energy-dependent term is the Laplace transform of their memory kernel
[Nakajima1958, Zwanzig1960](@cite). Projecting the assembled Liouvillian discards the factorization
of its dissipator into jump amplitudes, so a truncated effective Liouvillian obtained this way need
not be in GKSL form, exactly as for the Lie transform [Ikeda2021, Schnell2021](@cite). Chruściński
and Kossakowski showed that a Feshbach projection of the amplitudes, rather than of the map, yields
completely positive dynamics by construction [Chruscinski2013](@cite). See also
[CP-preserving completion](@ref cp-preserving-completion-theory).

## Comparison with the Lie transform

| | Lie transform | Bloch/Feshbach projection |
| --- | --- | --- |
| Solved for | micromotion ``\mathcal{K}`` | wave operator ``\Omega`` and ``\mathcal{B}`` |
| Normalization | ``\langle\mathcal{K}\rangle=0`` at every order | ``\langle\Omega\rangle=1``, then ``\Omega N=e^{\Lambda}`` with ``\langle\Lambda\rangle=0`` and ``\Lambda`` a series of nested commutators |
| Recurrence | nested commutators with lower-order ``\mathcal{K}`` | products of lower-order ``\Omega`` and ``\mathcal{B}`` |
| Result | van Vleck pair | the same van Vleck pair |

The projection trades the commutator recursion for a product recurrence followed by a
normalization. Usage is described in [Expansion algorithms](@ref expansion-algorithms-manual).

Mikami et al. use the Fourier convention of this page. Sources that write
``H(t)=\sum_m H_m e^{+im\omega t}``, such as Eckardt and Anisimovas, require ``m\mapsto-m``
[Eckardt2015](@cite).

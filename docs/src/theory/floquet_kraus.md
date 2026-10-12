```@meta
CurrentModule = FloquetExpansions
CollapsedDocStrings = true
```

# [Floquet–Kraus period map](@id floquet-kraus-theory)

The evolution of a periodically driven Lindblad system over one period is a completely positive,
trace-preserving (CPTP) map ``\Phi_T``. An effective generator describes it only on a branch where
``\Phi_T`` has a time-independent GKSL logarithm, and a truncated effective generator need not be
of GKSL form [Schnell2021, Ikeda2021](@cite). The period map itself has no such restriction. The
package therefore also expands ``\Phi_T`` directly and keeps it completely positive at every drive
frequency. As in the amplitude-space projection of [Chruscinski2013](@cite) and the jump-resolved
operator sum of [Reimer2019](@cite), complete positivity follows from truncating the amplitudes
before they are squared, and the exact map below is the Markovian limit of that operator sum. The
truncation is new: a joint expansion in the jump number and ``1/\omega_d`` in the van Vleck frame
of the no-jump evolution.

Throughout, ``T=2\pi/\omega_d``, the Fourier convention is ``A(t)=\sum_m A_m e^{-im\omega_d t}``,
and the order indexing is that of [High-frequency expansion](@ref high-frequency-expansion-theory).

## Quantum-jump amplitudes

For

```math
\dot\rho=-i[H(t),\rho]+\sum_a\gamma_a\mathcal D[L_a(t)]\rho,
```

the no-jump generator is ``Q(t)=-iH(t)-\tfrac12\sum_a\gamma_aL_a^\dagger(t)L_a(t)``, with
propagator ``U_0(t,s)=\mathcal T\exp\int_s^tQ``. Unravelling the dissipators gives the exact
period map as a sum over jump records,

```math
\Phi_T(\rho)=\sum_{q\ge0}\sum_{a_1\dots a_q}\gamma_{a_1}\cdots\gamma_{a_q}
\int_{0<t_1<\dots<t_q<T}K\rho K^\dagger\,dt_1\cdots dt_q,
```

```math
K=U_0(T,t_q)L_{a_q}(t_q)U_0(t_q,t_{q-1})\cdots L_{a_1}(t_1)U_0(t_1,0).
```

Each term is a positive square, so complete positivity is manifest. The amplitudes ``K`` form an
isometry: the map is trace preserving because the jump records exhaust every history.

## No-jump Floquet frame

The no-jump evolution is the Floquet problem of the non-Hermitian generator ``Q``. Its van Vleck
expansion gives

```math
U_0(t,s)=\mathcal P(\omega_dt)\,e^{(t-s)Q_{\rm eff}}\,\mathcal P(\omega_ds)^{-1},
```

with a periodic ``\mathcal P`` whose logarithm has zero period average. The Hori--Deprit and
Bloch/Feshbach algorithms return the same pair. A jump amplitude then reads

```math
K=\mathcal P(0)\,e^{(T-t_q)Q_{\rm eff}}\tilde L_{a_q}(t_q)\cdots
\tilde L_{a_1}(t_1)\,e^{t_1Q_{\rm eff}}\,\mathcal P(0)^{-1},
\qquad
\tilde L_a(t)=\mathcal P(\omega_dt)^{-1}L_a(t)\mathcal P(\omega_dt).
```

Between jumps the evolution is static, so only the transported jumps ``\tilde L_a`` carry
harmonics.

## Truncation and jump-time integrals

With ``t_j=Tu_j``, a record with ``q`` jumps contributes ``T^{q/2}`` to each amplitude. The static
exponentials are expanded too, since ``e^{T(u_{j+1}-u_j)Q_{\rm eff}}`` is a series in ``T``, so
every power of ``TQ_{\rm eff}`` and every further power of ``1/\omega_d`` in ``Q_{\rm eff}``,
``\mathcal P`` or ``\tilde L_a`` adds one order. Keeping the exponentials would resum all orders
and leave operator-valued jump-time integrals; expanding them makes every retained amplitude a
polynomial in ``1/\omega_d``. An amplitude of order ``q/2+n`` paired with one of order ``q/2+n'`` enters the map at
``\omega_d^{-(q+n+n')}``. Retaining every amplitude with ``q+n\le N`` therefore reproduces
``\Phi_T`` through ``\omega_d^{-N}``. The jump count is bounded only within one period; long-time
evolution still resums arbitrarily many jumps.

A retained amplitude is a finite sum ``\sum_\mu\mu(u)A_\mu`` of operators ``A_\mu`` times
monomials ``\mu(u)=\prod_ju_j^{\alpha_j}e^{-2\pi ik_ju_j}`` on the simplex
``0<u_1<\dots<u_q<1``. The jump-time integral pairs them through the Gram matrix

```math
G_{\mu\nu}=\int_{0<u_1<\dots<u_q<1}\mu(u)\overline{\nu(u)}\,du,
```

whose entries are exact rational functions of ``\pi``. Summing over jump records gives the
superoperator

```math
\Phi(\rho)=\sum_{q}\sum_{a_1\dots a_q}T^q\gamma_{a_1}\cdots\gamma_{a_q}
\sum_{\mu\nu}G_{\mu\nu}A_\mu\rho A_\nu^\dagger .
```

The output sideband of each jump stays an external label of the monomial. It is never projected
away, which is how the leading jump term reproduces ``\sum_mL_m\rho L_m^\dagger``.

## Complete positivity and trace preservation

The Choi matrix of each record is ``\int|K\rangle\langle K|\,du``, an integral of positive
operators, and ``G`` is its Gram form in the monomial basis. The retained amplitudes are squared
without truncation, so ``\Phi`` is completely positive for every drive frequency, not only
asymptotically. The amplitudes ``K(u)=\sum_\mu\mu(u)A_\mu`` are themselves a continuous Kraus
family of the map. Factoring ``G`` would give a finite Kraus form, but the exact pivots are
rational functions of ``\pi`` whose degree grows quickly with the order, so the package keeps the
division-free Gram form. The trace defect

```math
D=\Phi^\dagger(1)-1=\mathcal O(\omega_d^{-N-1})
```

vanishes through the retained order because the exact amplitudes form an isometry and every
omitted amplitude first enters at ``\omega_d^{-N-1}``. Beyond that order its sign is not fixed.
Every retained amplitude is itself truncated, so the map can increase the trace there, unlike the
truncation of [Reimer2019](@cite) in the jump number alone with exact no-jump propagators, which
is trace nonincreasing. The map ``\rho\mapsto\Phi(S\rho S)`` with
``S=(1+D)^{-1/2}`` is exactly trace preserving and still agrees with ``\Phi_T`` through
``\omega_d^{-N}``, wherever ``1+D`` is positive definite. This normalizes the amplitudes; it is
not a completion of a generator.

## Relation to effective generators

On a branch where the period map has a GKSL logarithm, the Floquet decomposition gives

```math
\Phi_T=\mathcal M(0)\,e^{T\mathcal L_{\rm eff}}\,\mathcal M(0)^{-1},
```

so ``T^{-1}\log\Phi_T`` is the stroboscopic generator at ``t=0``. It equals a van Vleck or GKSL
normal-form effective generator only after the micromotion ``\mathcal M(0)`` of that gauge is
undone; see [GKSL normal form](@ref gksl-normal-form-theory). A period map without a GKSL
logarithm is still a valid result.

Multi-jump records carry the physics that a direct Liouvillian expansion assigns to products of
dissipators. For a jump operator that rotates with the drive, such as
``L(t)=\sigma_z+\cos(\omega_dt)\sigma_y+\sin(\omega_dt)\sigma_x``, the van Vleck generator contains
the Hamiltonian shift ``-\tfrac{5\gamma^2}{4\omega_d}\sigma_z`` from the commutator of two
dissipator harmonics. In the period map the same term arises from ordinary two-jump amplitudes.

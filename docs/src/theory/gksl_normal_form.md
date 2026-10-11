```@meta
CurrentModule = FloquetExpansions
CollapsedDocStrings = true
```

# [GKSL normal form](@id gksl-normal-form-theory)

A truncated van Vleck expansion of a periodic Lindblad generator need not be of
Gorini--Kossakowski--Sudarshan--Lindblad (GKSL) form [Schnell2021, Ikeda2021](@cite). Positive
completion repairs a fixed truncation afterwards and changes only data beyond the retained order
[Haddadfarshi2015](@cite); see [CP-preserving completion](@ref cp-preserving-completion-theory).
The GKSL normal form instead uses a gauge freedom of the expansion itself. It chooses the static
part of the micromotion order by order so that every truncation of the effective generator is of
GKSL form by construction.

Throughout, ``\epsilon=\omega_d^{-1}``, ``[\mathcal A,\mathcal B]=\mathcal A\mathcal B-\mathcal
B\mathcal A`` is the commutator of linear maps, and the Fourier and order conventions are those of
[High-frequency expansion](@ref high-frequency-expansion-theory).

## Static gauge freedom

The homological equation fixes only the oscillatory part of each micromotion coefficient. The van
Vleck gauge sets the remaining static part to zero, ``\langle\mathcal K\rangle=0``. Keeping the
oscillatory part and adding a static similarity

```math
\mathcal W=e^{\mathcal S},
\qquad
\mathcal S=\sum_{n\ge1}\epsilon^n\mathcal S_n,
```

maps the effective generator to

```math
\mathcal L_{\rm eff}=\mathcal W^{-1}\mathcal L_{\rm VV}\mathcal W .
```

This is a formal similarity, so the spectrum agrees with the van Vleck spectrum order by order.
The coefficients do not. At order ``n`` the effective coefficient takes the form

```math
\mathcal E_n=\widehat{\mathcal V}_n+[\mathcal L_0,\mathcal S_n],
```

where ``\mathcal L_0`` is the period-averaged generator and ``\widehat{\mathcal V}_n`` collects
everything fixed by lower orders. Each ``\mathcal S_n`` is a free choice at its own order. It
enters the recurrence at all higher orders, so it is chosen inside the Bloch/Feshbach or
Hori--Deprit recurrence rather than applied to a finished expansion.

## Graded jump amplitudes

The effective generator is represented through channels born at definite orders. A channel ``a``
born at order ``m_a`` has rate ``\epsilon^{m_a}\gamma_a`` with ``\gamma_a>0`` and amplitude

```math
J_a(\epsilon)=l_a+\sum_{k\ge1}\epsilon^k\delta_a^{(k)} .
```

The input channels have ``m_a=0``. The finite generator at retained order ``N`` is

```math
\mathcal L^{[N]}
=-i\left[H^{[N]},\,\cdot\,\right]
+\sum_a\epsilon^{m_a}\gamma_a\,\mathcal D\!\left[J_a^{[N]}\right],
```

which is of GKSL form for every ``\epsilon`` since every rate is positive. In a dissipative frame
the Kossakowski coefficient of order ``n`` splits into three parts,

```math
c(\mathcal E_n)=K_n^<+K_n^{\rm lin}+P_n ,
```

with

```math
\begin{aligned}
K_n^<&=\sum_{a:\,m_a<n}\gamma_a\!\!\sum_{\substack{i+j=n-m_a\\ i,j\ge1}}\!\!
\delta_a^{(i)}\delta_a^{(j)\dagger},
\\
K_n^{\rm lin}&=\sum_{a:\,m_a<n}\gamma_a
\left(l_a\delta_a^{(n-m_a)\dagger}+\delta_a^{(n-m_a)}l_a^\dagger\right),
\\
P_n&=\sum_{a:\,m_a=n}\gamma_a\,l_a l_a^\dagger .
\end{aligned}
```

The first part is known from lower orders. The second is linear in the new corrections of the
active channels. The third holds the channels born at order ``n`` and is positive semidefinite.
Matching ``c(\mathcal E_n)`` at every retained order makes ``\mathcal L^{[N]}`` agree with the
expansion through ``\epsilon^N``; the untruncated products of the amplitudes supply a remainder
beyond that order.

## The static slot

The linear part reaches every Hermitian matrix whose block on the dark quotient vanishes, the
quotient of the frame by the span of the active amplitudes ``l_a``. Let ``Q`` project onto it
along that span. Only the dark block constrains the slot, and it gives the dark equation

```math
\Delta_n+\Phi_n(\mathcal S_n)=P_n\succeq0,
\qquad
\Delta_n=Q\left(c(\widehat{\mathcal V}_n)-K_n^<\right)Q^\dagger,
\qquad
\Phi_n(\mathcal S)=Q\,c\!\left([\mathcal L_0,\mathcal S]\right)Q^\dagger .
```

A solution ``\mathcal S_n`` fixes ``P_n``, whose factorization gives the newborn channels. The
remaining bright part of ``c(\mathcal E_n)-K_n^<-P_n`` then determines the corrections
``\delta_a^{(n-m_a)}``. Without dissipation there are no channels, ``\Delta_n`` vanishes, and
``\mathcal S_n=0`` returns the van Vleck coefficients.

The static gauge is searched among maps of GKSL shape,

```math
\mathcal S[\rho]=-i[h,\rho]
+\sum_{\mu\nu}s_{\mu\nu}\left(F_\mu\rho F_\nu^\dagger
-\tfrac12\{F_\nu^\dagger F_\mu,\rho\}\right),
```

with Hermitian ``h`` and ``s`` of either sign, built from frame operators ``F_\mu`` of bounded
polynomial degree. Families of increasing degree are tried in turn.

## Solving the dark equation

The dark equation is an affine positive-semidefinite feasibility problem. It is solved in exact
arithmetic in three steps:

1. **Facial reduction.** A principal block of ``\Delta_n`` that no gauge direction touches must
   already be positive semidefinite, and positivity then forces the rows over its kernel to
   vanish. These linear constraints are imposed until none remain. An indefinite untouched block
   proves that the family admits no solution.
2. **Canonical section.** Within the constrained affine family the least-squares point in the frame
   metric, of minimum gauge norm, is taken. If it is not positive semidefinite, rows of ``P_n`` are
   cancelled in frame order wherever the constraints stay consistent.
3. **Certificate.** An exact rate-weighted ``LDL^\dagger`` factorization certifies
   ``P_n=\sum_j\gamma_j l_j l_j^\dagger`` with rational rates and no square roots. Its columns are
   the newborn amplitudes.

The solution of the dark equation is not unique. The canonical section fixes one representative
given the frame order, and the spectrum through the retained order does not depend on the choice.

## Retained data and obstructions

The GKSL normal form and positive completion answer different questions. Completion keeps every
retained van Vleck coefficient and chooses only the omitted orders. The normal form keeps the
retained spectrum and changes the retained coefficients wherever some ``\mathcal S_n\neq0``,
compensating in the micromotion. For a driven Kerr oscillator with loss the normal form equals
the van Vleck coefficients through ``\epsilon^3``. At ``\epsilon^4`` a static gauge is needed and
a channel ``a^2-2a^\dagger a`` is born.

A static gauge of low polynomial degree need not exist at every order. For the same Kerr
oscillator the dark equation at ``\epsilon^6`` has no solution in any family of degree at most
three, so the normal form is available through ``\epsilon^5`` only. Such an obstruction is reported rather than approximated.

See [Floquet expansion](@ref) for the user workflow.

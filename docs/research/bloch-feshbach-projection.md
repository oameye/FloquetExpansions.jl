# Bloch/Feshbach projection for Floquet high-frequency expansions: literature research

Purpose: source material for rewriting `docs/src/theory/expansion_algorithms.md` into a brief theory page on the Bloch/Feshbach (wave-operator) route, parallel to `docs/src/theory/high_frequency_expansion.md`.

Evidence labels used throughout:

- **verified**: read in the primary source (section or equation given where possible).
- **abstract**: read only in the publisher abstract of the primary source.
- **secondary**: stated by another paper about the cited work; the citing paper is named.
- **inference**: our own reasoning, shown in place; **exact check** marks an inference also confirmed by exact rational arithmetic (Section 3.6).

## 1. Summary

1. The package recurrence is Bloch's energy-independent wave-operator equation on Sambe space with $P$ the zero Fourier sector. With $H_0=-M\omega$ that model space is exactly degenerate (eigenvalue 0), so Lindgren's generalized Bloch equation reduces to Bloch's degenerate form (inference; Lindgren form secondary via the Lindgren, Salomonson and Hedendahl 2009 preprint, `Lindgren2010`).
2. Mikami et al. (2016) already derive exactly this recurrence (their Eq. (22)). Their energy-independent "Brillouin–Wigner" Hamiltonian $H_{\rm BW}$ equals the package Bloch generator, $H_{\rm BW}=i\mathcal B$ (verified; identity of $B^{(2)}$ checked term by term).
3. Feshbach's $P/Q$ elimination gives the energy-dependent $H_{\rm eff}(\varepsilon)=PHP+PHQ(\varepsilon-QHQ)^{-1}QHP$ (Mikami Eq. (12), verified). The Bloch equation removes the energy dependence by replacing $\varepsilon\,\Omega$ with $\Omega\,H_{\rm B}$ (Mikami Eqs. (17) and (18), verified).
4. $H_{\rm B}$ is non-Hermitian because the projected eigenvectors are not orthonormal. The metric is $S=P\Omega^\dagger\Omega P=\langle\Omega^\dagger\Omega\rangle$ (Mikami Eq. (27); Kvaal 2008; Bravyi et al. 2011; verified).
5. des Cloizeaux: $H_{\rm dC}=S^{1/2}H_{\rm B}S^{-1/2}$. It equals the two-block canonical van Vleck and the Schrieffer–Wolff (direct-rotation) effective Hamiltonian (Kvaal Eq. (32) and Burgarth et al. Eq. (8.14), verified; Klein 1974, abstract).
6. **Key finding.** The Floquet van Vleck Hamiltonian is **not** $S^{1/2}H_{\rm B}S^{-1/2}$ in general. Exactly, $H_{\rm vV}=N^{-1}H_{\rm B}N$ with $N=\langle e^{-iK}\rangle=S^{-1/2}V$, $V$ unitary, so $H_{\rm vV}=V^\dagger H_{\rm dC}V$. Here $V=1+O(\omega^{-3})$. The two first differ at $\omega^{-3}$ when the harmonics allow $\langle (K^{(1)})^3\rangle\neq0$ (exact check). Jørgensen et al. (1975) and Mankodi and DiVincenzo (2024) report static analogues. Normalizations that remove different sets of couplings, or least action versus an off-diagonal generator, first diverge at third order in the perturbation.
7. The relation $F_{\rm vV}=N^{-1}H_{\rm BW}N$ with $N=\langle e^{-iK}\rangle$ follows directly from Mikami Eqs. (35), (38), (40) (inference from verified equations). The package conversion therefore has a Floquet precedent.
8. Liouvillians: every step is algebraic, and the similarity by $N$ is not unitary. By the argument in Section 6, item 6, TP and HP survive (inference), but CP does not (Burgarth et al. 2021 counterexample, verified; Kessler 2012, verified). Feshbach projection at the amplitude level (Chruściński and Kossakowski 2013) is the construction that keeps CPTP.
9. On cost, published statements are qualitative only. Bloch coefficients are quadratic polynomials in lower-order ones (Bravyi §1.2). Bloch's similarity transform scales linearly per order, while direct SW scales exponentially and optimized SW as $N^2$ (Araya Day et al. 2025, their classification). Mikami: Floquet–Magnus and van Vleck show "many more terms" than BW at high order (all verified). We found no published Floquet cost comparison.

## 2. Conventions map

Package: $\dot x=\mathcal G(t)x$, $\mathcal G(t)=\sum_m\mathcal G_m e^{-im\omega t}$, $\mathcal G=-iH$, micromotion $\mathcal M(t)=e^{\mathcal K(t)}$ with $\mathcal K=-iK$ for a Hermitian kick $K$. Sambe blocks $Q_{mn}=H_{m-n}-m\omega\delta_{mn}$ (Liouvillian: $i\mathcal L_{m-n}-m\omega\delta_{mn}$). Bloch objects: $\Omega(t)=\sum_m\Omega_m e^{-im\omega t}$, $\Omega_0=1$, $\partial_t\Omega=\mathcal G\Omega-\Omega\mathcal B$; Hamiltonian units $H_{\rm B}=i\mathcal B$. Order: $\omega^{-n}$ in component $n$.

| Source | Fourier / Sambe convention | Their objects | Translation to package |
| --- | --- | --- | --- |
| Mikami et al. 2016 (verified, Sec. II) | $H_{m,n}=\frac1T\int dt\,H(t)e^{i(m-n)\omega t}$, $T_m=\{e^{-im\omega t}\}$, $(H-M\omega)\lvert u\rangle=\varepsilon\lvert u\rangle$ | $P$ zero-photon ("photon vacuum"), $Q=1-P$, wave operator $\Omega(\varepsilon)$, $\Omega_{\rm BW}$, $H_{\rm eff}(\varepsilon)$, $H_{\rm BW}$, $\Xi(t)=\sum_m e^{-im\omega t}[\Omega_{\rm BW}]_{m,0}$, kick $\Lambda(t)$ with $U=e^{-i\Lambda(t)}e^{-iF(t-t')}e^{i\Lambda(t')}$ | Same sign as package: $H_{m,n}=H_{m-n}$. $H_{\rm BW}=H_{\rm B}=i\mathcal B$; $[\Omega_{\rm BW}]_{m,0}=\Omega_m$; $\Xi(t)=\Omega(t)$; $\Lambda=K$ |
| Eckardt and Anisimovas 2015 (verified, Secs. 2 to 4, App. C) | $H(t)=\sum_m e^{im\omega t}H_m$; $Q_{m'm}=H_{m'-m}+\delta_{m'm}m\omega$ | $\bar U_F=\exp\bar G$, $\bar G$ anti-Hermitian, block-off-diagonal in all photon sectors; $U_F(t)=e^{G(t)}$; kick $K=iG$; $H_F$ | $m\to-m$. $U_F(t)=\mathcal M(t)$, $G=\mathcal K$, $K_{\rm EA}=K$, $H_F=H_{\rm eff}$ |
| Bukov et al. 2015 (verified, Eqs. (1), (43)) | $H(t)=\sum_l H_le^{il\Omega t}$ | $U=e^{-iK(t_2)}e^{-iH_F(t_2-t_1)}e^{iK(t_1)}$ | $l\to-l$; kick sign as package |
| Rahav et al. 2003 (verified, Sec. III) | time domain; $\phi=e^{i\hat F}\psi$ | $\hat F$ periodic, static part "assumed to vanish" | $\hat F=K$ in the zero-average gauge |
| Bravyi et al. 2011 (verified, Secs. 1.2, 2, 3) | static, $H=H_0+\epsilon V$ | SW: $H_{\rm eff}=P_0UHU^\dagger P_0$, $U=e^{S}$; Bloch: $U=\begin{pmatrix}I&0\\U_{+-}&0\end{pmatrix}$, $H_{\rm eff}=P_0HU$ | Bloch $U$ = Sambe $\Omega$ with $P_0$ = zero sector |
| Kvaal 2008 (verified, Sec. III) | static | decoupling operator $\omega=Q\omega P$, $\Omega=P+\omega$; $H^{\rm BB}_{\rm eff}=PH(P+\omega)$; $H^{c}_{\rm eff}$ | $P+\omega^\dagger\omega=S$; $F=S^{1/2}$ |
| Lindgren, Salomonson, Hedendahl 2009 preprint of `Lindgren2010` (verified, Eqs. (1.1) to (1.10)) | static, $H=H_0+V$ | generalized Bloch eq. $[\Omega,H_0]P=(V\Omega-\Omega V_{\rm eff})_{\rm linked}P$, $V_{\rm eff}=PV\Omega P$ | Sambe: $H_0\to -M\omega$, $V\to$ Toeplitz $H$ |
| Viennot 2014 (verified, Sec. 2) | $i\hbar\,\partial_t$; initial-value problem | TDWO $\Omega(t)=U(t,0)(P_0U(t,0)P_0)^{-1}$, $i\hbar\dot\Omega=[H,\Omega]\Omega$, $H_{\rm eff}=P_0H\Omega$ | Not periodic; not the package $\Omega$ (see Q2) |
| Burgarth et al. 2021 (verified, Secs. I, VII to IX) | $e^{t(\gamma B+C)}$, arbitrary matrices | right and left wave operators $U_\ell$, $\tilde U_\ell$; $K_\ell=(\tilde U_\ell U_\ell)^{1/2}(\gamma B+D_\ell)(\tilde U_\ell U_\ell)^{-1/2}$ | $\tilde U U\to S$ in the unitary case |
| Keliri and Schirò 2026 (verified, Secs. II, III) | $\mathcal L(t)=\sum_m\mathcal L_m e^{-im\omega_0t}$; Sambe diagonal blocks $\mathcal L_0+im\omega_0$ | $L_{\rm eff}(\mu)$ by continued fractions | Same sign as package; their $L_F=-iQ$ (Liouvillian Sambe) |

## 3. Question 1: origins and definitions

### 3.1 Bloch wave operator and intermediate normalization

- The Bloch wave operator maps the model space onto the exact invariant subspace. It obeys $P\Omega=P$ (intermediate normalization) and $\Omega=\Omega P$. It is idempotent and solves $H\Omega=\Omega H\Omega$. The effective Hamiltonian is $H_{\rm B}=PH\Omega$, energy-independent, with the exact eigenvalues. Secondary via Viennot 2014, Eqs. (2) and (3) ("$\Omega^2=\Omega$ ... non-linear generalization of an eigenprojector"), and Bravyi et al. 2011 §1.2 ("$U$ is a (non-hermitian) projector onto the perturbed low-energy subspace", $[H_0+\epsilon V,U]U=0$). Both cite Bloch 1958. We could not read Bloch 1958.
- In Floquet form, the Sambe Bloch equation $Q\,\Omega=\Omega\,(PQ\Omega P)$ has $(m,0)$ component $\sum_nH_{m-n}\Omega_n-m\omega\Omega_m=\Omega_mH_{\rm B}$ with $H_{\rm B}=\sum_nH_{-n}\Omega_n=\langle H\Omega\rangle$. In the time domain this reads $i\partial_t\Omega=H\Omega-\Omega H_{\rm B}$, which is $\partial_t\Omega=\mathcal G\Omega-\Omega\mathcal B$ with $\mathcal B=-iH_{\rm B}$ (inference; matches Mikami Eq. (18) term by term).

### 3.2 Generalized (Lindgren) Bloch equation and linked diagrams

- Lindgren (1974) extends Rayleigh–Schrödinger perturbation theory to a model space "which is not necessarily degenerate" (abstract). The generalized Bloch equation $[\Omega,H_0]P=(V\Omega-\Omega V_{\rm eff})_{\rm linked}P$ with $V_{\rm eff}=PV\Omega P$ in intermediate normalization is secondary, via Lindgren, Salomonson and Hedendahl 2009, Eqs. (1.8) and (1.9). That chapter also reports the degenerate form $(E_0-H_0)\Omega P=(V\Omega-\Omega V_{\rm eff})_{\rm linked}P$ attributed to Bloch, Eq. (1.1).
- For Floquet: $H_0=-M\omega$ and $P$ = zero sector give $[\Omega,H_0]P$ with component $m\omega\Omega_m$. This is the degenerate case $E_0=0$ (inference).
- Brandow (1967) rederives the Goldstone expansion "starting from Brillouin–Wigner (BW) perturbation theory" (abstract). Via Lindgren et al. 2009 (secondary): Brandow extends the linked-diagram theorem to quasi-degeneracy with folded diagrams.
- Lindgren et al. 2009 (secondary, for Lindgren 1978 and Mukherjee): for a complete model space, the normal-ordered exponential $\Omega=\{\exp T\}$ has a completely connected $T$.

### 3.3 Why $H_{\rm B}$ is non-Hermitian

- Kvaal 2008 §III.C (verified): $H^{\rm BB}_{\rm eff}=PH(P+\omega)$ is non-Hermitian because "rejecting the excluded space eigenvector components renders the effective eigenvectors non-orthonormal", $\langle\psi_j^{\rm eff}|\psi_k^{\rm eff}\rangle=\delta_{jk}-\langle\psi_j|Q|\psi_k\rangle$.
- Mikami 2016 Eq. (27) (verified), Floquet: the zero-photon projections are not orthonormal, and orthonormality reads $u^{0\dagger}_\alpha[\Omega_{\rm BW}^\dagger\Omega_{\rm BW}]_{0,0}u^0_\beta=\delta_{\alpha\beta}$. Mikami also shows that at $n$th order the non-Hermitian part changes eigenvalues and right eigenvectors only at $O(\omega^{-n-1})$ (Sec. I).
- Package: $iB$ is non-Hermitian from $\omega^{-2}$ (exact check, Section 3.6).

### 3.4 des Cloizeaux and the canonical (Hermitian) effective Hamiltonian

- Burgarth et al. 2021 (verified, Sec. I; secondary for des Cloizeaux 1960): "des Cloizeaux showed that one can turn the non-skew-Hermitian $\gamma B+D$ into a skew-Hermitian $\gamma B+K$ by an additional similarity transformation keeping the block structure."
- Kvaal 2008 Eqs. (22), (31), (32) (verified): $F=(P+\omega^\dagger\omega)^{1/2}$ and $H^{c}_{\rm eff}=FH^{\rm BB}_{\rm eff}F^{-1}$. Here $H^{c}_{\rm eff}=Pe^{-G}He^{G}P$ with $e^G$ the direct rotation of $P$ onto the exact subspace. Kvaal attributes $G=\tanh^{-1}(\omega-\omega^\dagger)$ to Shavitt and Redmon 1980 (secondary).
- Klein 1974 (abstract): the formulations "due to Van Vleck, Kemble, and Primas, to des Cloizeaux, and to Buleavski" are shown to yield identical results.
- Burgarth Eqs. (8.13), (8.14) (verified): in the unitary case the generalized SW map reduces to the Bravyi et al. direct rotation, and $K_\ell=(\tilde U_\ell U_\ell)^{1/2}(\gamma B+D_\ell)(\tilde U_\ell U_\ell)^{-1/2}$ with $\tilde U_\ell=U_\ell^\dagger$.
- Package translation: $S=P\Omega^\dagger\Omega P=\sum_m\Omega_m^\dagger\Omega_m=\langle\Omega^\dagger\Omega\rangle$ (Parseval), and $H_{\rm dC}=S^{1/2}H_{\rm B}S^{-1/2}$.

### 3.5 Feshbach projection, Rayleigh–Schrödinger vs Brillouin–Wigner, and removal of the energy dependence

- Feshbach P/Q elimination gives the energy-dependent effective Hamiltonian $H_{\rm eff}(z)=PHP+PHQ(z-QHQ)^{-1}QHP$. Verified as a formula in Durand and Paidarová 2006, Eq. (6), and in Mikami Eq. (12) (Floquet form $H_{\rm eff}(\varepsilon)=PH[1-\tfrac{Q}{\varepsilon+M\omega}H]^{-1}P$). Attribution to Feshbach 1958 and 1962 is secondary via Chruściński and Kossakowski 2013, ref. [20]; we could not read Feshbach's papers.
- Bravyi §1.2 (verified): this self-energy form yields a $z$-dependent $H_{\rm eff}(z)$, and "there is no systematic way of getting rid of this dependence".
- Mikami Sec. II.B (verified) gives the mechanism of removal. Using $PH\Omega(\varepsilon)P=H_{\rm eff}$, one replaces $\Omega(\varepsilon)\varepsilon$ on the right by $\Omega(\varepsilon)PH\Omega(\varepsilon)$ in Eq. (17). This yields the $\varepsilon$-independent Eq. (18), $\Omega=P+\frac{Q}{M\omega}H\Omega-\frac{Q}{M\omega}\Omega PH\Omega$, and its $1/\omega$ recursion Eq. (22). The BW denominators $1/(\varepsilon+m\omega)$ become RS-type $1/(m\omega)$ with the effective Hamiltonian acting from the right. Mikami also expands the BW series in $1/\omega$ to recover the same $H_{\rm BW}$ (text after Eq. (23f)).
- Naming: Mikami call the $\varepsilon$-independent object $H_{\rm BW}$. It is Bloch's energy-independent effective Hamiltonian in intermediate normalization. The theory page should say so, to avoid readers expecting a self-consistent BW calculation.
- Durand and Paidarová 1998 (abstract): energy-dependent wave-operator theory for time-independent and time-dependent Hamiltonians, relation to the $(t,t')$ method, three-space partition. Durand 1983 (abstract): the wave operator solves an equation "which generalizes those previously established by Bloch, Löwdin, Jørgensen, and Lindgren".

### 3.6 Exact check of the normalizations (inference, exact rational arithmetic)

Setup: random $3\times3$ Hermitian $H_0$, random $H_{\pm m}$ with $H_{-m}=H_m^\dagger$, entries in $\{-3,\dots,3\}+i\{-3,\dots,3\}$, `Complex{Rational{BigInt}}`, orders $\omega^0$ to $\omega^{-6}$. The script implements the package Bloch recurrence, $S=\langle\Omega^\dagger\Omega\rangle$, $H_{\rm dC}=S^{1/2}H_{\rm B}S^{-1/2}$ by binomial series, and $N$ fixed order by order from $\langle\log(\Omega N)\rangle=0$, with $H_{\rm vV}=N^{-1}H_{\rm B}N$. Results:

- $\Lambda=\log(\Omega N)$ is anti-Hermitian (so $e^\Lambda$ is unitary) with $\langle\Lambda\rangle=0$; $H_{\rm vV}$ and $H_{\rm dC}$ are Hermitian at all orders, and $H_{\rm B}$ is non-Hermitian from $\omega^{-2}$.
- The harmonics $k\neq0$ of $\Omega(t)^\dagger\Omega(t)$ vanish, so $\Omega(t)^\dagger\Omega(t)=S$ for every $t$. $NN^\dagger=S^{-1}$ exactly, and $V=S^{1/2}N$ is unitary.
- The first order at which $V\neq1$, and $H_{\rm vV}\neq H_{\rm dC}$, depends on the harmonics: $\omega^{-3}$ for $\{\pm1,\pm2\}$, $\omega^{-4}$ for $\{\pm1,\pm3\}$, and $\omega^{-5}$ for $\{\pm1\}$ only.

The mechanism: with $\langle K\rangle=0$, $N=\langle e^{-iK}\rangle=1-\tfrac12\langle K^2\rangle+\tfrac i6\langle K^3\rangle+\dots$. The first anti-Hermitian term is $\tfrac{i}{6}\langle (K^{(1)})^3\rangle\omega^{-3}$. It needs three harmonics summing to zero, so it vanishes when only odd harmonics are present, which is consistent with the table above.

## 4. Question 2: Floquet and time-dependent formulations

| Work | P/Q onto zero-photon sector? | Formulation and names |
| --- | --- | --- |
| Mikami et al. 2016 | **Yes** (verified). "Project the whole Hilbert space $\mathbb H\otimes\mathbb T$ onto the zero-photon subspace"; $[P]_{mn}=\delta_{mn}\delta_{m0}$ | "wave operator" $\Omega(\varepsilon)$, $\Omega_{\rm BW}$; $H_{\rm eff}(\varepsilon)$ (BW, self-consistent); $H_{\rm BW}$ (energy-independent); $\Xi(t)$ labelled "wave op." in Fig. 1, "a counterpart to the so-called micromotion operator" |
| Eckardt and Anisimovas 2015 | **No two-block P/Q.** Multi-partition with projectors $\bar P_m$ on every photon sector (App. C, Eqs. (C.5) to (C.9)) | "canonical van Vleck degenerate perturbation theory" (citing Shavitt and Redmon 1980); $\bar U_F=\exp\bar G$, $\bar G$ block-off-diagonal "In order to minimize the mixing"; uniqueness from $\bar G_D=0$ (C.19); $K=iG$ "kick operator" |
| Rahav, Gilary, Fishman 2003 | No. Time-domain order-by-order gauge transformation | $\hat F$ periodic, static part set to zero |
| Bukov, D'Alessio, Polkovnikov 2015 | No. Time domain, van Vleck HFE vs Floquet–Magnus; §4.2 compares with SW | "kick operator" $K(t)$; SW "coincides with the HFE at least up to the second order, and agrees with the ME up to a static gauge transformation" (verified) |
| Jolicard and Killingbeck 2003 (Part II) | Active space $P_0$, general. Extended-space wave-operator equation (abstract) | Two equations of motion: "a non-linear differential equation in the usual Hilbert space" and one "in an extended Hilbert space with an extra time variable ... equivalent to the usual Bloch equation when the Floquet Hamiltonian is taken in place of the ordinary Hamiltonian" (abstract) |
| Viennot 2014 | Active space $P_0$, initial-value problem | TDWO $\Omega(t)=P(t)(P_0P(t)P_0)^{-1}$, $\Omega(0)=P_0$; Eq. (8) $(H-i\hbar\partial_t)\Omega=\Omega(H-i\hbar\partial_t)\Omega$, "a Bloch equation with the Floquet Hamiltonian ... in the extended Hilbert space" (verified) |
| Jolicard et al. 1994 (JCP 100, 325) | Active space of Floquet (quasivibrational) states for photodissociation | Bloch wave operator "to select the active space" (abstract). Not a $1/\omega$ expansion |
| Sadreev 2012 | No zero-photon projection. $P$ = driven wire ⊗ Floquet index, $Q$ = leads | Floquet index as extra lattice dimension; "non-Hermitian effective Hamiltonian" (abstract) |
| Keliri and Schirò 2026 (preprint) | **Yes**, Liouville–Sambe (verified, Sec. III.A) | $L_{\rm eff}(\mu)=\mathcal L_0+\mathcal L_-T_1^+(\mu)\mathcal L_++\mathcal L_+T_{-1}^-(\mu)\mathcal L_-$ by matrix continued fractions: an energy-dependent Feshbach/BW form, self-consistent |

Notes:

- The Jolicard/Viennot time-dependent wave operator is an initial-value object, $\Omega(0)=P_0$ with $P(t)$ the evolved projector. It is not the periodic, zero-average-normalized $\Omega(t)$ of the package. The package $\Omega$ is the Sambe Bloch operator of the extended-space equation that the Jolicard–Killingbeck abstract names (inference).
- Mikami Eq. (39) gives the time-domain equation $-i\partial_t\Xi^{-1}=\Xi^{-1}H-H_{\rm BW}\Xi^{-1}$. It is equivalent to the package Bloch equation for $\Xi=\Omega$ (inference).
- Mikami Eqs. (35), (38), (40) (verified) are $F_{\rm FM}=e^{-i\Lambda_{\rm vV}(t_0)}F_{\rm vV}e^{i\Lambda_{\rm vV}(t_0)}$, $F_{\rm FM}=\Xi(t_0)H_{\rm BW}\Xi(t_0)^{-1}$, and $\Xi(t)^{-1}=\frac1T\int_0^Tdt'\,e^{-i\Lambda(t')}e^{i\Lambda(t)}$. Hence $\Xi(t)=e^{-i\Lambda(t)}\langle e^{-i\Lambda}\rangle^{-1}$, and $F_{\rm vV}=\langle e^{-i\Lambda}\rangle^{-1}H_{\rm BW}\langle e^{-i\Lambda}\rangle$. This is the package $\mathcal G_{\rm eff}=N^{-1}\mathcal BN$, $N=\langle e^{\mathcal K}\rangle$ (inference from verified equations).

## 5. Question 3: comparison with the other methods

### 5.1 Schrieffer–Wolff

- SW 1966 (verified): a canonical transformation $e^{S}He^{-S}$ with $S$ first order in $V$, chosen so that $[H_0,S]=H_1$ eliminates $V_{kd}$ "to first order".
- Bravyi et al. 2011 (verified): SW is the direct rotation $U=\sqrt{R_{P_0}R_P}$ (Def. 2.2). Lemma 2.3 gives uniqueness: there is a unique anti-Hermitian $S$ with $e^{S}Pe^{-S}=P_0$, $S$ block-off-diagonal ($P_0SP_0=Q_0SQ_0=0$), and $\|S\|<\pi/2$. The perturbative series (Eqs. (3.4) to (3.13)) is built from nested commutators with $\tanh(x/2)$ and $x\coth x$ Taylor coefficients.
- The relation to Bloch: $H_{\rm SW}=H_{\rm dC}=S^{1/2}H_{\rm B}S^{-1/2}$ in the two-block case (Kvaal; Burgarth, verified).

### 5.2 Canonical van Vleck / quasidegenerate perturbation theory

- Shavitt and Redmon 1980 (abstract): three quasidegenerate theories are compared "in terms of a common general formulation based on a similarity transformation which decouples the model space and complementary space components". Their canonical van Vleck form is compared "with the approach based on intermediate normalization". Content beyond the abstract is secondary via Kvaal.
- Jørgensen, Pedersen, Chedin 1975 (verified, §§2, 7). Effective Hamiltonians from different choices of which couplings to remove are related by a unitary within the model space (Eq. (2.13)). The anti-Hermitian part of $p=PUP$ "determines $U$ and $A$ completely" (§7, Theorem). The choice $\{p\}_A=0$ gives Soliverez's operators, and different choices agree at orders 1 and 2 (Eq. (7.11)). "For the near-degenerate case deviations occur in third order."
- Mankodi and DiVincenzo 2024 (preprint, abstract): for multi-block diagonalization, the "least action" condition and block-off-diagonality of the generator "diverge at third order".
- Araya Day et al. 2025 (verified, §§1, 3.2, 3.4): multi-block SW generalizations "yield different effective Hamiltonians when applied to more than two subspaces". In the $2\times2$ case, $U_{AA}=U_{AA}^\dagger$, $U_{BB}=U_{BB}^\dagger$, $U_{AB}=-U_{BA}^\dagger$ is equivalent to SW.

### 5.3 Sambe-space block diagonalization (Floquet)

- Eckardt and Anisimovas 2015 (verified): the unitaries are translationally invariant in the photon index, and correspond to time-periodic $U(t)$. The generator is block-off-diagonal with respect to all photon sectors, which in the time domain is $\langle G\rangle=0$. The result is the package van Vleck pair (the high_frequency_expansion.md claim).
- Inference: the des Cloizeaux/SW rotation for the two-block split (zero sector vs rest) is not Toeplitz. Its generator has no $QQ$ block, while a Toeplitz generator with nonzero $(0,n)$ blocks also has nonzero $(1,n+1)$ blocks. Its restriction to $P$ is still the zero column of a periodic unitary $\mathcal M_{\rm dC}(t)=\Omega(t)S^{-1/2}=e^{-iK(t)}V^\dagger$. This is a valid Floquet micromotion in another static gauge: $\langle\mathcal M_{\rm dC}\rangle=S^{-1/2}$ is Hermitian positive rather than $\langle\log\mathcal M\rangle=0$. This matches the Jørgensen/Soliverez condition $\{p\}_A=0$.

### 5.4 Lie transform / Hori–Deprit

- Hori 1966 (verified, §1): a canonical map through a Lie series $f(x,y)=\sum_n\frac{\varepsilon^n}{n!}D_S^nf$, with $D_S$ the Poisson bracket with $S$.
- Deprit 1969 (verified, abstract and §1): Lie transforms with an $\varepsilon$-dependent generator, recursive substitution through "explicit chains of Poisson brackets".
- Package: $\mathcal K$ solved directly through the homological equation, with $\langle\mathcal K\rangle=0$. The result coincides with the Sambe multi-block canonical van Vleck, and with $\Lambda=\log(\Omega N)$ from the Bloch route (package fact checked by the caller; consistent with Section 3.6).

### 5.5 Hamiltonian vs non-Hermitian statements

- Unitary-only: des Cloizeaux ($S^{1/2}$ needs an adjoint), SW direct rotation, Hermiticity of $H_{\rm vV}$ and $H_{\rm dC}$, and $\Omega(t)^\dagger\Omega(t)=S$.
- Carries over algebraically: the Bloch equation and recurrence (Burgarth: arbitrary matrices, verified). Also the van Vleck condition $\langle\Lambda\rangle=0$, the similarity $N^{-1}\mathcal BN$, and uniqueness of the zero-average normal form (package).
- The non-Hermitian analogue of $S$ is $\tilde UU$ built from left and right wave operators (Burgarth Eqs. (7.2), (8.2), (8.14)). The package does not need it, because its canonical gauge is $\langle\Lambda\rangle=0$, not a metric condition.

## 6. Question 4: open systems

1. Chruściński and Kossakowski 2013 (verified): Feshbach projection at the amplitude level on $\mathcal H_S\otimes\mathcal H_E$ with a pure environment reference state. Eq. (13), $\Lambda_t(\lvert\psi\rangle\langle\psi\rvert)=Z_t\lvert\psi\rangle\langle\psi\rvert Z_t^\dagger+{\rm Tr}_E(Y_t\cdot Y_t^\dagger)$, yields a CPTP map. Their Born-like approximation "leads to the legitimate completely positive and trace preserving quantum evolution", which the Nakajima–Zwanzig route does not guarantee. The method is not periodic-drive specific.
2. Burgarth et al. 2021 (verified): Bloch's theory generalized to arbitrary generators, with a symmetric similarity that reduces to des Cloizeaux in the unitary case. $D$, $\tilde D$, $K$ are HP and TP (Sec. IX); CP "is not guaranteed". A counterexample (Sec. X C) shows an ideal effective generator (same block structure, similar, near identity, CP) "cannot always be provided".
3. Kessler 2012 (verified): SW for Liouvillians is a "non-unitary rotation". The second order is of Lindblad form; the third order is HP and TP but "not obviously of Lindblad form". Using the gauge freedom to enforce Lindblad form is left to future work.
4. Keliri and Schirò 2026 (preprint, verified): the zero-photon Feshbach $L_{\rm eff}(\mu)$ is TP and HP. They add "there is no guarantee that $L_F$ will be of Lindblad form" (Sec. II).
5. Sadreev 2012 (abstract): Feshbach with the Floquet index gives a non-Hermitian effective Hamiltonian for open leads (a scattering problem, not a Lindbladian).
6. Package-level inference (Liouvillian Bloch): if every harmonic is TP ($\langle\!\langle 1\rvert\mathcal L_m=0$), then $\langle\!\langle1\rvert\mathcal R^{(n)}_0=0$ because $\Omega^{(j)}_0=0$ for $j\ge1$, so every $\mathcal B^{(n)}$ is TP. With conjugation $J(A)=A^\dagger$ and $J\mathcal L_mJ=\mathcal L_{-m}$, induction on the recurrence gives $J\Omega^{(n)}_mJ=\Omega^{(n)}_{-m}$ and $J\mathcal B^{(n)}J=\mathcal B^{(n)}$, i.e. HP. $\Lambda$ is a Lie series in the $\mathcal L_m$, so $\langle\!\langle1\rvert\Lambda=0$ and $\langle\!\langle1\rvert N=\langle\!\langle1\rvert$, and $\mathcal G_{\rm eff}$ is TP and HP. CP is not preserved (items 2 to 4; Schnell 2021 and Ikeda 2021 already in `refs.bib`). Not checked numerically.

## 7. Question 5: computational aspects (only what sources say)

- Bravyi et al. 2011 §1.2 (verified): "The advantage of the Bloch expansion is that the perturbative series for $H_{\rm eff}$ has somewhat simpler structure compared with the corresponding expansion in the SW method ... one can express $U_n$ as a second-degree polynomial in $U_1,\dots,U_{n-1}$."
- Araya Day et al. 2025 §§1, 3.2 (verified): SW recursions "require an exponentially growing number of matrix products in each order". "The direct computation of the series elements requires $\sim\exp N$ multiplications, and even an optimized one has a $\sim N^2$ scaling." They list Bloch's similarity transform among the algorithms meeting their requirement of Cauchy-product ($\sim N$ per order) scaling without products by $H_0$, which "preserves the Hamiltonian spectrum while breaking its hermiticity".
- Mikami et al. 2016 §II.D (verified): "there appear many more terms (especially if the commutators are expanded) at higher order in the Floquet–Magnus and van Vleck expansions than in the BW expansion ... higher-order terms can be efficiently computed in the BW expansion by using the simple recursion relation (22)." This is qualitative, with no operation counts.
- Nikolaev 2016 §V (verified): for static canonical PT to order 32, the Kato-based explicit algorithm is "faster ... than the Magnus expansion, but less efficient then the Van Vleck method". "The comparison of efficiency ... is not unambiguous." This concerns Kato, not Bloch.
- Connectedness:
  - Bravyi Theorem 2 (verified) is a linked-cluster theorem in the spatial sense for the SW $H_{\rm eff}$ on lattices. Via additivity, it "applies to any operator-valued multivariate series that obeys the additivity property".
  - Lindgren (secondary via Lindgren et al. 2009): linked-diagram form of the generalized Bloch equation and a connected $T$ in $\{\exp T\}$ for complete model spaces.
  - Inference for the package: $\Lambda$ is a Lie series in the harmonics (it equals the Hori–Deprit $\mathcal K$). The recurrence terms $\Omega^{(j)}\mathcal B^{(n-j)}$ are model-space returns, the analogue of Brandow's folded diagrams.
  - We found no source that states the free-Lie (Lyndon) structure of $\log\Omega$ for Floquet Bloch operators.
- Package-internal (not literature): issues #96 and #69 count $(N-1)(N+2)/2$ series products for the Bloch core through package order $N$, explicitly not a runtime claim.

## 8. Comparison table

| | Bloch (intermediate) | Feshbach / BW | des Cloizeaux / SW (two-block) | Canonical van Vleck, Sambe multi-block (EA 2015) | Hori–Deprit (package) |
| --- | --- | --- | --- | --- | --- |
| Solved for | wave operator $\Omega=P+X$, $X=Q\Omega P$ | resolvent $(\varepsilon-QHQ)^{-1}$, self-consistent $\varepsilon$ | unitary $e^S$ ($S$ off-diagonal P/Q) or $\Omega$ plus $S^{\pm1/2}$ | Toeplitz unitary $e^{\bar G}$ | periodic $\mathcal K(t)$ |
| Normalization / uniqueness | $P\Omega P=P$ ($\langle\Omega\rangle=1$) | none (energy-dependent) | $S^\dagger=-S$, $PSP=QSQ=0$, $\lVert S\rVert<\pi/2$ (Bravyi L2.3) | $\bar G^\dagger=-\bar G$, $\bar G_D=0$ in all photon sectors (EA C.19) | $\langle\mathcal K\rangle=0$ |
| Hermiticity of $H_{\rm eff}$ (Hamiltonian) | non-Hermitian from $\omega^{-2}$ | Hermitian for real $\varepsilon$, but $\varepsilon$-dependent | Hermitian | Hermitian | Hermitian |
| Operations per order | products, quadratic in lower orders | products and resolvents | nested commutators (tanh/coth weights) | nested commutators, extra multi-block terms (EA C.13) | nested commutators (BCH triangles) |
| Effective operator | $H_{\rm B}=\langle H\Omega\rangle=i\mathcal B$ | $H(\varepsilon)$; $H_{\rm B}$ after removing $\varepsilon$ | $H_{\rm dC}=S^{1/2}H_{\rm B}S^{-1/2}$ | $H_{\rm vV}=N^{-1}H_{\rm B}N=V^\dagger H_{\rm dC}V$ | $=H_{\rm vV}$ |
| Periodic micromotion | $\Omega(t)$, not unitary; $\Omega^\dagger\Omega=S$ | none directly | $\Omega(t)S^{-1/2}=e^{-iK}V^\dagger$ (unitary) | $e^{-iK(t)}$ | $e^{-iK(t)}$ |
| Liouvillian version | yes (Burgarth; package) | yes (Keliri; CK at amplitude level) | needs a left wave operator (Burgarth) | similarity, not unitary | similarity, not unitary |

Relations, with $N=\langle e^{\mathcal K}\rangle=\langle e^{-iK}\rangle=S^{-1/2}V$: $H_{\rm vV}=N^{-1}H_{\rm B}N$ and $H_{\rm dC}=S^{1/2}H_{\rm B}S^{-1/2}$. $V=1+O(\omega^{-3})$, with the first nonzero order drive dependent (Section 3.6).

## 9. Suggested outline for the theory page ("Bloch/Feshbach projection")

Mirrors `high_frequency_expansion.md`: one intro paragraph, six sections, and an optional comparison. It should describe no API and link to the manual page `expansion-algorithms-manual` for usage.

0. **Intro.**
   - The fast Fourier sectors can be eliminated by projection instead of transformation. The result is the same van Vleck pair after a static normalization.
   - Cite [Bloch1958, Feshbach1962, Killingbeck2003, Jolicard2003, Mikami2016Brillouin].
1. **Model space in Sambe space.**
   - Define $P$ = zero Fourier sector, $Q=1-P$; relate to the $Q_{mn}$ blocks of `floquet_theory.md`.
   - With $H_0=-M\omega$, $P$ is an exactly degenerate model space.
   - Cite [Sambe1973, Mikami2016Brillouin, Lindgren1974].
2. **Feshbach elimination and the energy dependence.**
   - One equation: $H_{\rm eff}(\varepsilon)=H_0+\sum_{m,n\neq0}H_{-m}[(\varepsilon+M\omega-QHQ)^{-1}]_{mn}H_n$.
   - One sentence: its eigenvalues are self-consistent solutions.
   - Cite [Feshbach1958, Feshbach1962, Mikami2016Brillouin]. Optionally [Brandow1967] for BW-based linked expansions.
3. **Bloch equation and recurrence.**
   - The wave operator with $P\Omega P=P$ and the time-domain equation $\partial_t\Omega=\mathcal G\Omega-\Omega\mathcal B$, with $\mathcal B=\langle\mathcal G\Omega\rangle$.
   - The recurrence and $\Omega^{(n+1)}_m=\frac im\mathcal R^{(n)}_m$.
   - State the mechanism: $\varepsilon$ is replaced by right multiplication with $\mathcal B$.
   - Say that Mikami's $H_{\rm BW}$ is this $i\mathcal B$.
   - Cite [Bloch1958, Lindgren1974, Mikami2016Brillouin, Bravyi2011].
4. **Leading terms and non-Hermiticity.**
   - $\mathcal B^{(0)}$, $\mathcal B^{(1)}$ (equal to van Vleck) and $\mathcal B^{(2)}$ in the package form (issue #96).
   - The metric $S=\langle\Omega^\dagger\Omega\rangle$ and why $H_{\rm B}$ is non-Hermitian.
   - des Cloizeaux's $S^{1/2}H_{\rm B}S^{-1/2}$, noting its equality with two-block SW and canonical van Vleck.
   - Cite [Bloch1958, DesCloizeaux1960, Kvaal2008, Klein1974, ShavittRedmon1980, Bravyi2011].
5. **Normalization to the van Vleck gauge.**
   - $\Omega N=e^\Lambda$, $\langle\Lambda\rangle=0$, $N=\langle e^\Lambda\rangle$, $\mathcal G_{\rm eff}=N^{-1}\mathcal BN$, $N^{(2)}$.
   - For Hamiltonians, $N=S^{-1/2}V$, so the van Vleck and des Cloizeaux representatives differ by a model-space unitary. They first differ at $\omega^{-3}$ or later depending on the harmonics.
   - $\Lambda$ is a Lie series; one sentence on the connected (commutator) structure.
   - Cite [Mikami2016Brillouin] for the Floquet relation, [Jorgensen1975] for the unitary freedom within the model space, [Eckardt2015, horiTheory1966, Deprit1969].
6. **Liouvillians.**
   - Same equations with $\mathcal G=\mathcal L$. Similarity, not unitarity. TP and HP retained, CP not.
   - Amplitude-level Feshbach as the CP-preserving counterpart. Link to CP-preserving completion.
   - Cite [Burgarth2021, Kessler2012, Chruscinski2013, Schnell2021, Ikeda2021]; optionally [Keliri2026].
7. **Comparison** (optional short table or two sentences).
   - Products vs commutators; cite [Bravyi2011, ArayaDay2025].
   - Keep qualitative, with no complexity claims beyond the sources.

## 10. BibTeX entries not yet in `docs/src/refs.bib`

Checked against the current `refs.bib` keys; none of these keys exist. Metadata from Crossref records (DOI resolver) or the arXiv listing.

```bibtex
@article{Kvaal2008,
  author = {Kvaal, Simen},
  title = {Geometry of effective {Hamiltonians}},
  journal = {Physical Review C},
  volume = {78},
  number = {4},
  pages = {044330},
  year = {2008},
  doi = {10.1103/PhysRevC.78.044330},
  eprint = {0808.1831},
}

@article{Bravyi2011,
  author = {Bravyi, Sergey and DiVincenzo, David P. and Loss, Daniel},
  title = {{Schrieffer}--{Wolff} transformation for quantum many-body systems},
  journal = {Annals of Physics},
  volume = {326},
  number = {10},
  pages = {2793--2826},
  year = {2011},
  doi = {10.1016/j.aop.2011.06.004},
}

@article{ShavittRedmon1980,
  author = {Shavitt, Isaiah and Redmon, Lynn T.},
  title = {Quasidegenerate perturbation theories. {A} canonical van {Vleck} formalism and its relationship to other approaches},
  journal = {The Journal of Chemical Physics},
  volume = {73},
  number = {11},
  pages = {5711--5717},
  year = {1980},
  doi = {10.1063/1.440050},
}

@article{Klein1974,
  author = {Klein, D. J.},
  title = {Degenerate perturbation theory},
  journal = {The Journal of Chemical Physics},
  volume = {61},
  number = {3},
  pages = {786--798},
  year = {1974},
  doi = {10.1063/1.1682018},
}

@article{Lindgren1974,
  author = {Lindgren, I.},
  title = {The {Rayleigh}-{Schrödinger} perturbation and the linked-diagram theorem for a multi-configurational model space},
  journal = {Journal of Physics B: Atomic and Molecular Physics},
  volume = {7},
  number = {18},
  pages = {2441--2470},
  year = {1974},
  doi = {10.1088/0022-3700/7/18/010},
}

@article{Brandow1967,
  author = {Brandow, Baird H.},
  title = {Linked-Cluster Expansions for the Nuclear Many-Body Problem},
  journal = {Reviews of Modern Physics},
  volume = {39},
  number = {4},
  pages = {771--828},
  year = {1967},
  doi = {10.1103/RevModPhys.39.771},
}

@article{Jorgensen1975,
  author = {Jørgensen, Flemming and Pedersen, Thorvald and Chedin, Alain},
  title = {A projector formulation for the {Van Vleck} transformation. {III}. {Generalization} and relation to the contact transformation},
  journal = {Molecular Physics},
  volume = {30},
  number = {5},
  pages = {1377--1395},
  year = {1975},
  doi = {10.1080/00268977500102911},
}

@article{Burgarth2021,
  author = {Burgarth, Daniel and Facchi, Paolo and Nakazato, Hiromichi and Pascazio, Saverio and Yuasa, Kazuya},
  title = {Eternal adiabaticity in quantum evolution},
  journal = {Physical Review A},
  volume = {103},
  number = {3},
  pages = {032214},
  year = {2021},
  doi = {10.1103/PhysRevA.103.032214},
  eprint = {2011.04713},
}

@article{Chruscinski2013,
  author = {Chruściński, Dariusz and Kossakowski, Andrzej},
  title = {{Feshbach} Projection Formalism for Open Quantum Systems},
  journal = {Physical Review Letters},
  volume = {111},
  number = {5},
  pages = {050402},
  year = {2013},
  doi = {10.1103/PhysRevLett.111.050402},
}

@article{Sadreev2012,
  author = {Sadreev, Almas F.},
  title = {{Feshbach} projection formalism for transmission through a time-periodic potential},
  journal = {Physical Review E},
  volume = {86},
  number = {5},
  pages = {056211},
  year = {2012},
  doi = {10.1103/PhysRevE.86.056211},
}

@article{Kessler2012,
  author = {Kessler, E. M.},
  title = {Generalized {Schrieffer}-{Wolff} formalism for dissipative systems},
  journal = {Physical Review A},
  volume = {86},
  number = {1},
  pages = {012126},
  year = {2012},
  doi = {10.1103/PhysRevA.86.012126},
  eprint = {1205.5440},
}

@article{ArayaDay2025,
  author = {Araya Day, Isidora and Miles, Sebastian and Kerstens, Hugo and Varjas, Daniel and Akhmerov, Anton R.},
  title = {{Pymablock}: An algorithm and a package for quasi-degenerate perturbation theory},
  journal = {SciPost Physics Codebases},
  pages = {50},
  year = {2025},
  doi = {10.21468/SciPostPhysCodeb.50},
  eprint = {2404.03728},
}

@article{Durand1983,
  author = {Durand, Philippe},
  title = {Direct determination of effective {Hamiltonians} by wave-operator methods. {I}. {General} formalism},
  journal = {Physical Review A},
  volume = {28},
  number = {6},
  pages = {3184--3192},
  year = {1983},
  doi = {10.1103/PhysRevA.28.3184},
}

@article{Durand1998,
  author = {Durand, Philippe and Paidarová, Ivana},
  title = {Wave operator theory of quantum dynamics},
  journal = {Physical Review A},
  volume = {58},
  number = {3},
  pages = {1867--1878},
  year = {1998},
  doi = {10.1103/PhysRevA.58.1867},
}

@article{Viennot2014,
  author = {Viennot, David},
  title = {Almost quantum adiabatic dynamics and generalized time-dependent wave operators},
  journal = {Journal of Physics A: Mathematical and Theoretical},
  volume = {47},
  number = {6},
  pages = {065302},
  year = {2014},
  doi = {10.1088/1751-8113/47/6/065302},
  eprint = {1308.1528},
}

@article{Nikolaev2016,
  author = {Nikolaev, Andrey},
  title = {{Kato} expansion in quantum canonical perturbation theory},
  journal = {Journal of Mathematical Physics},
  volume = {57},
  number = {6},
  pages = {062102},
  year = {2016},
  doi = {10.1063/1.4953639},
  eprint = {1504.05113},
}

@misc{Keliri2026,
  author = {Keliri, Andriani and Schirò, Marco},
  title = {{Sambe} Approach to {Floquet}-{Lindblad} Open Quantum Systems},
  year = {2026},
  eprint = {2606.09727},
  archiveprefix = {arXiv},
  doi = {10.48550/arXiv.2606.09727},
}

@misc{Mankodi2024,
  author = {Mankodi, Ishan N. H. and DiVincenzo, David P.},
  title = {Perturbative power series for block diagonalisation of {Hermitian} matrices},
  year = {2024},
  eprint = {2408.14637},
  archiveprefix = {arXiv},
  doi = {10.48550/arXiv.2408.14637},
}
```

Optional (metadata verified via Crossref; content **not** read):

```bibtex
@article{Okubo1954,
  author = {Ôkubo, Susumu},
  title = {Diagonalization of {Hamiltonian} and {Tamm}-{Dancoff} Equation},
  journal = {Progress of Theoretical Physics},
  volume = {12},
  number = {5},
  pages = {603--622},
  year = {1954},
  doi = {10.1143/PTP.12.603},
}

@article{Kato1949,
  author = {Kato, T.},
  title = {On the Convergence of the Perturbation Method. {I}},
  journal = {Progress of Theoretical Physics},
  volume = {4},
  number = {4},
  pages = {514--523},
  year = {1949},
  doi = {10.1143/ptp/4.4.514},
}

@article{Jolicard1994,
  author = {Jolicard, Georges and Killingbeck, John P. and Durand, Philippe and Heully, Jean Louis},
  title = {A wave operator description of molecular photodissociation processes using the {Floquet} formalism},
  journal = {The Journal of Chemical Physics},
  volume = {100},
  number = {1},
  pages = {325--333},
  year = {1994},
  doi = {10.1063/1.467001},
}
```

Secondary sources used above as verification (cite only if the page needs them):

```bibtex
@article{Durand2006,
  author = {Durand, Philippe and Paidarová, Ivana},
  title = {From effective {Hamiltonians} to fluctuation and dissipation},
  journal = {Theoretical Chemistry Accounts},
  volume = {116},
  number = {4-5},
  pages = {559--565},
  year = {2006},
  doi = {10.1007/s00214-006-0101-9},
}

@incollection{Lindgren2010,
  author = {Lindgren, Ingvar and Salomonson, Sten and Hedendahl, Daniel},
  title = {Coupled Clusters and Quantum Electrodynamics},
  booktitle = {Recent Progress in Coupled Cluster Methods},
  series = {Challenges and Advances in Computational Chemistry and Physics},
  publisher = {Springer},
  pages = {357--374},
  year = {2010},
  doi = {10.1007/978-90-481-2885-3_13},
}
```

Caveats on the BibTeX:

- `Lindgren2010`: the text read was the authors' preprint dated 2009-05-12 (fy.chalmers.se/~f3ail/Publications/CCQED.pdf). Equation numbers quoted above as "Lindgren et al. 2009" refer to that preprint. The published pages and DOI come from Crossref, and the editors were not checked.
- Crossref gives Kato's volume 4 issue 4 with year 1949; Semantic Scholar lists 1951.

## 11. Open questions and unverified items

1. **Not read in full (claims are secondary):**
   - Bloch 1958, des Cloizeaux 1960, Feshbach 1958 and 1962: publisher paywall; no OA copy per Unpaywall.
   - Shavitt and Redmon 1980; Klein 1974: abstract only.
   - Killingbeck and Jolicard 2003 Parts I and II: abstracts only. Bot protection blocked the PDF even though Unpaywall lists it as OA; the Chrome extension was not connected.
   - Brandow 1967, Lindgren 1974, Durand 1983 and 1998, Sadreev 2012: abstracts only.
   - Okubo 1954, Kato 1949: nothing beyond metadata.
2. The des Cloizeaux formula $S^{1/2}H_{\rm B}S^{-1/2}$ is verified in Kvaal (Eqs. (31), (32)) and Burgarth (Eq. (8.14)), not in des Cloizeaux 1960 itself.
3. Mikami et al. has an erratum, PRB 99, 019902 (2019), DOI 10.1103/PhysRevB.99.019902. Its content was not accessible; arXiv v3 (May 2016) predates it. Check whether it affects Eqs. (22), (23), (38) to (40) before citing equation numbers.
4. Crossref lists a one-page item, J. Phys. A 39, 13591 (2006), DOI 10.1088/0305-4470/39/43/C01, under the Part I title. It is probably a corrigendum to Killingbeck and Jolicard 2003 Part I; not checked.
5. We found no Floquet source that states $H_{\rm vV}\neq S^{1/2}H_{\rm B}S^{-1/2}$. The claim rests on our exact check plus the static analogues (Jørgensen 1975; Mankodi and DiVincenzo 2024; Araya Day 2025). Order counting differs: their "third order" is in the perturbation strength, while $\omega^{-3}$ is fourth order in $H$. Issue #97's wording "up to an allowed unitary rotation inside the model space" is consistent. The caller's example "$H_{\rm vv}=S^{1/2}H_BS^{-1/2}$" is not, as an exact identity.
6. `high_frequency_expansion.md`, "Sambe-space block diagonalization", says "Seek a block-off-diagonal generator $S$". For the package van Vleck gauge this must mean block-off-diagonal in all photon sectors (EA C.19) with Toeplitz structure. Off-diagonality in $P$ vs $Q$ alone gives des Cloizeaux, which can differ from order $\omega^{-3}$ on (drive dependent). This is worth one clarifying phrase there, but it is outside this page's scope.
7. The Liouvillian TP/HP statements for $\mathcal B$, $N$, $\mathcal G_{\rm eff}$ (Section 6, item 6) are inference and were not checked numerically.
8. We found no published operation count comparing the Floquet Bloch recurrence with the Hori–Deprit recursion. The package counts (#96, #69) are internal.

## 12. Points in the current draft `expansion_algorithms.md` that the sources contradict or refine

1. The draft says the split into the zero sector and its complement "is Feshbach's projection formalism". The P/Q split is shared by every method. Feshbach's result is the energy-dependent $H(\varepsilon)$. The implemented recurrence is Bloch's energy-independent equation (Lindgren's form, degenerate case). The page should state both and give the mechanism, $\varepsilon\to$ right multiplication by $\mathcal B$ (Mikami Eqs. (17), (18)).
2. The draft says "The same projection underlies the Brillouin–Wigner high-frequency expansion". This understates it. Mikami's Eq. (22) is the package recurrence, $H_{\rm BW}=i\mathcal B$, and Mikami's Eqs. (35), (38), (40) already contain $F_{\rm vV}=N^{-1}H_{\rm BW}N$ with $N=\langle e^{-iK}\rangle$. The normalization section should cite Mikami for the relation.
3. The draft cites only [Bloch1958, DesCloizeaux1960] for non-Hermiticity. Add Mikami Eq. (27) for the Floquet case. If des Cloizeaux is mentioned as the Hermitian fix, the page must add that in general it is not the van Vleck gauge from $\omega^{-3}$ on.
4. The draft table row "Normalization: none" for Hori–Deprit could mislead. Hori–Deprit imposes $\langle\mathcal K\rangle=0$ directly; the Bloch route imposes $\langle\Omega\rangle=1$ and converts afterwards.
5. The remaining draft statements check out: $\mathcal B=\langle\mathcal G\Omega\rangle$, the recurrence, $\mathcal B^{(1)}$, $N^{(2)}$, $\mathcal G^{(2)}_{\rm eff}=\mathcal B^{(2)}+[\mathcal G_0,N^{(2)}]$, the Lie-series claim, and no CP guarantee. Each is consistent with the sources and the exact check.

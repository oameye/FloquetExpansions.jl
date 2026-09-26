> Notice: this file summarizes two unpublished research notes of the maintainer; it must not be committed or published without the maintainer's decision.

# Floquet Nakajima–Zwanzig and Feshbach/Bloch notes: review for the `BlochFeshbach()` theory page

## Sources and scope

- **NZ note**: `Open-quantum-system-Floquet-expansion/Floquet-Nakajima-Zwanzig/floquet_nakajima_zwanzig.tex` (665 lines). Sections used most: §1 "The Nakajima–Zwanzig projection" (lines 105–156), §2 "Born–Markov reduction and the vertex expansion" (158–239), §4 "The slow invariant subspace and all-order equivalence" (251–326), §5 TCL (328–354). §§6–16 (positivity, self-consistent and constrained schemes, CP expansion, Pawula boundary) are read for classification only.
- **CP note**: `.../cp_feshbach_vanvleck_full_note.tex` (1775 lines). Sections used most: "The Floquet NZ reduction already contains Van Vleck" (83–134), "Brillouin–Wigner: Floquet Feshbach projection as a high-frequency expansion" (136–221), the open-leg Bloch–Feshbach/Hori–Deprit recurrences (615–766), "Feshbach/Bloch and canonical Van Vleck are gauge choices of the same slow subspace" (930–965), and "FloquetExpansions.jl certification record" (1446–1518).
- Package context read: `docs/src/theory/{high_frequency_expansion,floquet_theory,expansion_algorithms,cp_completion}.md`, `docs/agents/documentation.md`, `src/{bloch_projection,bloch_connected_van_vleck,expansion_algorithms,engine}.jl`, `test/expansion_algorithms.jl`, `CONTEXT.md`, issues #69, #72, #97, #102, #168, PR #122 and PR #323.

Line numbers below refer to the `.tex` files; equation labels are the notes' own `\label`s.

## 1. Convention dictionary

The notes and the package differ in the Fourier sign and in order indexing. Every translation below was checked against the package recurrence numerically (Section 5).

| Item | Notes | Package | Translation |
| --- | --- | --- | --- |
| Fourier sign | $\mathcal L(t)=\sum_m\mathcal L_m e^{+im\omega t}$ (NZ eq:nz-harmonics, line 109; CP line 88; CP §§ 967–1227 implicitly) | $\mathcal G(t)=\sum_m\mathcal G_m e^{-im\omega t}$ | $\mathcal L^{\rm notes}_m=\mathcal G_{-m}$ |
| Exceptions | CP lines 636, 682, 787, 859, 868, 1478 use $e^{-im\omega t}$ ("FE.jl's convention") | | no relabel there |
| Sambe operator | $\mathbb L=\mathcal L(\theta)-\omega\partial_\theta$, acting as $\mathcal L_0-im\omega$ on harmonic $m$ of $\varrho(\theta)=\sum\varrho_m e^{im\theta}$ (NZ eq:nz-sambe, eq:nz-split) | $\mathbb G_{mn}=\mathcal G_{m-n}+im\omega\,\delta_{mn}=-iQ_{mn}$, with $Q_{mn}=H_{m-n}-m\omega\delta_{mn}$ (`floquet_theory.md:107,119`) | $\mathbb L=\mathbb G$ after $m\to-m$; energy-like $Q=i\mathbb G$ |
| Energy-like BW operator | $\mathbb K=H(\theta)-i\omega\partial_\theta$ (CP line 140) | $Q$ | $\mathbb K=Q$ after relabel; CP eq:BWresolvent, $(\epsilon-H_0+m\omega)^{-1}$, already matches the package sign $Q_{mm}=H_0-m\omega$ |
| Graph / wave operator | $P+Y$, $PYP=0$ (NZ eq:nz-invariance) | periodic $\Omega(t)$, $\langle\Omega\rangle=1$; Sambe column $(\Omega_m)$ | $Y_m^{\rm notes}=\Omega_{-m}$ |
| Slow generator | $\mathcal L_{\rm eff}$ (NZ/Bloch) | Bloch effective generator $\mathcal B$ (internal, `bloch.effective`) | $\mathcal L_{\rm eff}=\mathcal B$ |
| van Vleck generator | $\mathcal L_{\rm VV}$ | $\mathcal G_{\rm eff}$ (`effective_component`) | $\mathcal L_{\rm VV}=\mathcal G_{\rm eff}$ |
| Normalization | $N=Pe^{G}P$, $\mathcal L_{\rm eff}=N\mathcal L_{\rm VV}N^{-1}$ (NZ eq:nz-gauge-N; CP eq:NZVV) | $N=\langle e^{\Lambda}\rangle$, $\mathcal G_{\rm eff}=N^{-1}\mathcal B N$ | identical statement |
| Micromotion generator | Sambe $G$ with $PGP=0$, $G_m\simeq\mathcal L_m/(im\omega)$ | $\mathcal K=\Lambda$, $\mathcal K^{(1)}_m=(i/m)\mathcal G_m$ | $G^{\rm notes}_m=\mathcal K_{-m}$ |
| Order index (NZ note) | $\mathcal L^{(1)}_{\rm eff}=\mathcal L_0$, $\mathcal L^{(2)}\propto\omega^{-1}$, $\mathcal L^{(3)}\propto\omega^{-2}$ (lines 199–207) | component $n$ carries $\omega^{-n}$; `order` $N$ keeps $n=0,\dots,N-1$ (`src/engine.jl:26-27`) | notes' order $k$ = package component $k-1$; package `order` $N$ = notes' orders $1..N$ |
| Order index (CP note) | $\mathcal L^{(1)}_{\rm VV}\propto\omega^{-1}$, $\mathcal L^{(2)}_{\rm VV}\propto\omega^{-2}$ (eq:vv1, eq:vv2), but $\mathcal L_{\rm CP}^{[N+1]}$ in NZ counting (eq:LCP) | as above | mixed within the CP note |
| Open-leg grading | $\tau=\omega t$, $A_\delta=\delta A_1+\delta^2A_2$, $\delta=\omega^{-1/2}$ (CP eq:openleg-A) | $\omega^{-n}$ grading in $t$ units | for a single grading $A=\mathcal G/\omega$: $B^{\rm open}_{n+1}=\mathcal B^{(n)}$, $X_n=\Omega^{(n)}$ |
| Vertex count | $k$ = number of nonzero-harmonic insertions, with $\mathcal L_0$ kept exact inside resolvents (NZ §2, §4) | not used; $\mathcal G_0$ enters the residual perturbatively | the package series is the re-expansion of the vertex series in $\operatorname{ad}_{\mathcal G_0}/(m\omega)$ |
| Units | $\hbar=1$, frequency units | same | none |

## 2. Question 1: Bloch/Feshbach content of the notes, in package conventions

### 2.1 Sambe lift and P/Q partition

NZ §1 (eqs nz-sambe to nz-slow/fast) and CP §2 lift $\dot\rho=\mathcal L(t)\rho$ to the autonomous Sambe generator and split it with the period-average projector $P\varrho=\varrho_0$, $Q=1-P$. In package form: $\mathbb G=\mathbb G_0+\mathbb G_c$ with $\mathbb G_0$ harmonic-diagonal ($\mathcal G_0+im\omega$ on sector $m$) and $\mathbb G_c$ the nonzero harmonics; $P\mathbb G_0P=\mathcal G_0$, $P\mathbb G_cP=0$.

### 2.2 Invariance (graph) equation = Bloch equation

NZ eq:nz-invariance and eq:nz-sylvester (lines 261–268), CP eq:graphL: $\mathbb L(P+Y)=(P+Y)\mathcal L_{\rm eff}$, $\mathcal L_{\rm eff}=\mathcal L_0+P\mathbb L_cY$. In package form this is exactly the Sambe form of the Bloch equation $\partial_t\Omega=\mathcal G\Omega-\Omega\mathcal B$:

$$\mathbb G\,(P+X)=(P+X)\,\mathcal B,\qquad X=(\Omega_m)_{m\ne0},$$

- P row: $\mathcal B=\mathcal G_0+\sum_{n\ne0}\mathcal G_{-n}\Omega_n=\langle\mathcal G\Omega\rangle$.
- Q row, after eliminating $\mathcal B$ (NZ eq:nz-riccati, CP eq:Riccati), a Riccati equation:

$$(\operatorname{ad}_{\mathcal G_0}+im\omega)\,\Omega_m=-\mathcal G_m-\sum_{n\ne0,m}\mathcal G_{m-n}\Omega_n+\Omega_m\,(\mathcal B-\mathcal G_0),\qquad m\ne0 .$$

Consequences stated in NZ §4 (lines 278): the graph of $X$ is an invariant subspace of $\mathbb G$; every eigenvalue of $\mathcal B$ is an exact Floquet exponent (the branch that connects to $X=0$); $\Omega_m$ reconstructs the micromotion pointwise. The notes identify the intermediate normalization $P\Omega P=P$ with the time-averaging projector (CP line 935, $P\Omega_BP=P$), i.e. $\langle\Omega\rangle=1$.

### 2.3 Homological hierarchy, vertex expansion, Born–Markov

NZ eq:nz-homological (lines 285–290) expands the Riccati equation in the vertex count $k$ with $\mathcal L_0$ exact. Package form:

$$(\operatorname{ad}_{\mathcal G_0}+im\omega)X^{(k)}_m=-\mathcal G_m\delta_{k1}-\sum_{n\ne0,m}\mathcal G_{m-n}X^{(k-1)}_n+\sum_{j=1}^{k-2}X^{(j)}_m\sum_{n\ne0}\mathcal G_{-n}X^{(k-1-j)}_n .$$

The two-vertex truncation is the Born–Markov generator (NZ eq:nz-bm-ad), in package form $\mathcal B_{\rm BM}=\mathcal G_0-\sum_{m\ne0}\mathcal G_{-m}(\operatorname{ad}_{\mathcal G_0}+im\omega)^{-1}\mathcal G_m$. Expanding $-(\operatorname{ad}+im\omega)^{-1}=\frac{i}{m\omega}-\frac{\operatorname{ad}}{(m\omega)^2}+\dots$ gives the package $\mathcal B^{(1)}$ and the two-vertex part of $\mathcal B^{(2)}$ (NZ lines 199–203). The three-vertex term is NZ eq:nz-three-vertex; its leading part is the three-harmonic double sum.

The package recurrence is this hierarchy regraded by $\omega^{-n}$ with $\mathcal G_0$ moved into the residual. Explicitly, from the recurrence (checked numerically):

$$\mathcal B^{(2)}=-\sum_{m\ne0}\frac{\mathcal G_{-m}[\mathcal G_0,\mathcal G_m]}{m^2}-\sum_{m\ne0}\sum_{n\ne0,m}\frac{\mathcal G_{-m}\mathcal G_{m-n}\mathcal G_n}{mn}.$$

### 2.4 Brillouin–Wigner and the energy-independent wave operator

CP §3 (eq:BWFeshbach, eq:BWwave, eq:BWresolvent, eq:LiouvillianBW) restates Mikami et al.: $H_{\rm BW}(\epsilon)=P\mathbb KP+P\mathbb KQ(\epsilon-Q\mathbb KQ)^{-1}Q\mathbb KP$ on the zero-photon sector, intermediate normalization $P\Omega_{\rm BW}P=P$, and expansion of $(\epsilon-H_0+m\omega)^{-1}$ in $1/(m\omega)$ as the source of the high-frequency series. The same form is written for a Floquet Liouvillian, $\mathcal L_{\rm BW}(z)=P\mathbb LP+P\mathbb LQ(z-Q\mathbb LQ)^{-1}Q\mathbb LP$. The note's key sentence (lines 210–219): the graph $Y$ is the energy-independent, operator-valued counterpart of the BW wave operator; the Riccati equation constructs the whole slow subspace at once instead of one root $z$ at a time. Summary chain eq:FeshbachBWVVchain: spectral Feshbach ↔ BW, invariant graph ↔ Bloch wave operator, canonical normalization ↔ van Vleck (see Section 6, item F5, for a caveat on the last link).

Package form of the root relation: if $\mathcal Bv=\lambda v$ then $\mathcal B_{\rm BW}(\lambda)v=\lambda v$ with $\mathcal B_{\rm BW}(z)=P\mathbb GP+P\mathbb GQ(z-Q\mathbb GQ)^{-1}Q\mathbb GP$ (checked numerically).

### 2.5 Normalization to van Vleck

NZ eq:nz-vv-graph, eq:nz-gauge-N, eq:nz-gauge-N-expanded (lines 310–324); CP eq:NZVV (131), eq:BlochVV (958–963). Uniqueness of the perturbative graph gives $e^{G}P=(P+Y)N$ with $N=Pe^GP$ and $\mathcal L_{\rm eff}=N\mathcal L_{\rm VV}N^{-1}$. Expansion: $N=1+\sum_{m\ne0}\mathcal L_m\mathcal L_{-m}/(2(m\omega)^2)+O(\omega^{-3})$, invariant under $m\to-m$, hence identical to the draft's $N^{(2)}=\frac12\sum_{m\ne0}\mathcal G_m\mathcal G_{-m}/m^2$ (`expansion_algorithms.md:135`). This is the gauge that reconciles the NZ third order with the van Vleck third order (NZ line 203). For a Hamiltonian, CP eq:BWVVsimilarity (lines 188–200) introduces the metric $S=P\Omega^\dagger\Omega P=\langle\Omega^\dagger\Omega\rangle$ and $H=S^{1/2}H_{\rm BW}S^{-1/2}$ "up to a unitary rotation inside the model space".

### 2.6 Open-leg Bloch–Feshbach recurrence

CP eqs openleg-wave to openleg-homological (lines 654–686): $A_\delta\Omega-\partial_\tau\Omega=\Omega B$, $\mathcal R_n=A_1X_{n-1}+A_2X_{n-2}-\sum_{j=1}^{n-1}X_jB_{n-j}$, $B_n=P\mathcal R_n$, $(X_n)_\ell=(i/\ell)(\mathcal R_n)_\ell$ in the $e^{-i\ell\tau}$ convention. With a single grading $A=\mathcal G/\omega$ this is literally the package recurrence (`src/bloch_projection.jl:155-177`), with the index shift of Section 1. CP eq:openleg-B22 gives $-i\sum_{\ell>0}[A_{1,\ell},A_{1,-\ell}]/\ell$, which is the package $\mathcal B^{(1)}=i\sum_{m>0}[\mathcal G_{-m},\mathcal G_m]/m$. CP eq:openleg-HD-residual has the same Lie weights as the engine (`src/engine.jl:67-76,115-127`: $(-1)^j/j!$ and $(-1)^{j}/(j+1)!$ with a minus sign). The open-leg content itself (output fields, Kraus grades) is later-stack material.

### 2.7 Trace sum rule and spectral exactness

NZ lines 232–237: $\operatorname{Tr}\mathcal L_{\rm eff}=\operatorname{Tr}\mathcal L_0$ for every scheme; for commutator series because commutators are traceless, for the ad-resolvent generator by the pairing identity eq:nz-trace-sum-rule. For the package this means $\operatorname{Tr}\mathcal B^{(n)}=\operatorname{Tr}\mathcal G_{\rm eff}^{(n)}=0$ for $n\ge1$ (similarity preserves the formal trace). Checked numerically.

## 3. Question 2: Nakajima–Zwanzig as Feshbach, and TCL

- **What the notes say.** NZ §1 applies Nakajima–Zwanzig to the Sambe-lifted equation with $P$ the period average (the "irrelevant" part is micromotion, not a bath). The exact memory equation is NZ eq:nz-exact with kernel $\mathcal K(\tau)=P\mathbb L_ce^{Q\mathbb LQ\tau}\mathbb L_cP$ (eq:nz-memory-kernel, also CP eq:NZkernel). Its Laplace transform is the Feshbach self-energy $P\mathbb L_c(z-Q\mathbb LQ)^{-1}\mathbb L_cP$, i.e. the P/Q Feshbach partition of CP eq:LiouvillianBW in the time domain (the notes do not write this Laplace step explicitly; it is standard). The Markov closure $\bar\sigma(t-\tau)=e^{-\mathcal L_{\rm eff}\tau}\bar\sigma(t)$ (eq:nz-markov) is shown in NZ §4 to be the Sylvester pair eq:nz-sylvester, i.e. the Bloch invariance equation. So: NZ kernel (time domain) ↔ Feshbach/BW self-energy (Laplace domain, energy dependent) ↔ Bloch wave operator (energy independent, algebraic).
- **TCL.** NZ §5: with the undressed period-average projector the exact TCL generator $\mathcal K_{\rm TCL}(t)=\dot G(t)G(t)^{-1}$, $G(t)=Pe^{\mathbb Lt}P$, is exact but quasi-periodic and can become singular; projecting instead onto $\operatorname{ran}(P+Y)$ (the Bloch subspace, oblique spectral projector) collapses it to the constant $\mathcal L_{\rm eff}$ (lines 354). TCL at second order averages to the Born–Markov generator (eq:nz-tcl2).
- **Recommendation for the page.** At most two sentences in the Liouvillian paragraph: with $P$ the period average, the P/Q partition of the Sambe-lifted master equation is the Nakajima–Zwanzig projection; its memory kernel is the time-domain Feshbach self-energy, and the Bloch generator is the constant generator on the invariant subspace. Cite Nakajima 1958, Zwanzig 1960 and a textbook (Breuer–Petruccione). Leave out the Born–Markov, self-consistent, TCL, and memory-kernel extensions: they are not implemented, and the "exact Markov closure" wording rests on integral representations that do not converge as written (item F3).

## 4. Question 3: consistency with the implementation and the draft

### 4.1 Agreements (all checked, Section 5)

1. The notes' invariance equation, translated, is the package Bloch equation; the notes' $\mathcal L_{\rm eff}$ is the package's internal $\mathcal B$ and $\langle\Omega\rangle=1$ is the notes' $P\Omega P=P$.
2. The package recurrence (`src/bloch_projection.jl:155-177`, draft lines 77–85) is the $\omega^{-n}$ expansion of the notes' graph solution: truncation error $O(\omega^{-(N+1)})$ against the exact invariant subspace.
3. $\mathcal B^{(1)}$ equals the van Vleck first order (draft line 96; NZ eq:nz-second-order-fourier after $m\to-m$; CP eq:VV1-package).
4. $N=\langle e^\Lambda\rangle$, $\mathcal K=\Lambda$, $\mathcal G_{\rm eff}=N^{-1}\mathcal BN$ (draft lines 106–130; `src/expansion_algorithms.jl:88-108`) agree with NZ eq:nz-gauge-N and CP eq:BlochVV. The draft's $N^{(2)}$ and $\mathcal G^{(2)}_{\rm eff}=\mathcal B^{(2)}+[\mathcal G_0,N^{(2)}]$ (draft lines 135–138) agree with NZ eq:nz-gauge-N-expanded and reproduce the notes' third-order formula eq:nz-third-order-fourier = CP eq:vv2, which is invariant under $m\to-m$.
5. Hamiltonian bookkeeping: inverse weight $1/m$ for $H$ and $i/m$ for $\mathcal L$ (`src/expansion_algorithms.jl:4-5`), log phase $1$ vs $i$ and kick phase $i$ vs $1$ (lines 185–186, 203–204) are consistent with $\mathcal G=-iH$, $\mathcal K=-iK$; CP eq:BlochHom's leading $Y_m=H_m/(m\omega)$ matches $1/m$ once its sign is set to the package convention (item in Section 1).
6. Non-Hermiticity of the Bloch Hamiltonian (draft line 101; CP lines 186, 1399) is confirmed: $H_B-H_B^\dagger=O(\omega^{-2})$.
7. The draft's statement that the implementation expands $\Lambda$ in nested commutators with rational coefficients (lines 144–149) matches `bloch_compile_connected_log_words` and `bloch_lyndon_decomposition` (`src/bloch_connected_van_vleck.jl:146-267,327-353`).
8. The rotating-jump coefficient $-5\gamma^2/(4\omega)\,\sigma_z$ (CP eq:RR-fivefourths) holds in the package convention.

### 4.2 Contradictions and hazards

- **C1 (draft vs notes, hard).** `docs/src/theory/expansion_algorithms.md:13-14` says the Bloch/Feshbach route "solves a linear equation for a wave operator", and lines 87–89 say "The recurrence is linear". The notes (NZ eq:nz-riccati, line 282; CP eq:Riccati, line 120) show the wave-operator equation is quadratic (Riccati) once $\mathcal B=\langle\mathcal G\Omega\rangle$ is inserted; the fold term $\sum_j\Omega^{(j)}\mathcal B^{(n-j)}$ (`src/bloch_projection.jl:162-164`) is that quadratic back-action. What is true: each order is a linear triangular step (as it is for Hori–Deprit), and the step uses ordinary products instead of nested commutators.
- **C2 (hazard for the rewrite).** CP eq:BWVVsimilarity and eq:FeshbachBWVVchain ("canonical normalization ↔ Van Vleck") must not be read as "the package's van Vleck generator is $S^{1/2}H_BS^{-1/2}$". Numerically the package's zero-average representative differs from the isometric (des Cloizeaux) one by a static unitary: $\|H_{\rm VV}-S^{1/2}H_BS^{-1/2}\|\cdot\omega^3\approx0.84$ at $\omega=120,240,480$ for a qubit drive with harmonics $\pm1,\pm2$ (and about $\omega^{-5}$ for a single harmonic). Correspondingly $N=\langle e^\Lambda\rangle$ is not Hermitian ($N-N^\dagger=O(\omega^{-3})$), so $N\ne S^{-1/2}$. The draft as written does not claim otherwise.
- **C3 (notes vs code, stale).** CP line 968 ("the simpler construction currently implemented in the package: the coherent-frame/one-dissipator reduction") and CP lines 398, 781: PR #122 is closed unmerged (base `research/104-cp-hfe-oracle`); no transport code exists on `main` or on this branch. See Section 5b.
- **C4 (notes, convention hazard).** CP eq:vv1 (line 1038), $\sum_{m>0}\frac{i}{m\omega}[\mathcal L_m,\mathcal L_{-m}]$, is in the $e^{+im\omega t}$ convention, while CP eq:VV1-package (line 879), $-\frac{i}{\omega}\sum_{m>0}[\mathcal L_m,\mathcal L_{-m}]/m$, is in the package convention. Both are right; the note does not restate the convention between them.
- **C5 (notes vs code, scope of a cost claim).** CP line 448: the Bloch–Feshbach wave-operator recurrence has a quadratic order-table product count. True for the recurrence (series level $(N-1)(N+2)/2$, issue #69), but the implemented `BlochFeshbach()` also runs the connected-log conversion, whose product count grows about fourfold per order for harmonics $\{0,\pm1,\pm2\}$: log-bracket products 16, 84, 358, 1482, 6228, 26796 for `order` 3 to 8, against 60 to 560 harmonic-level recurrence generator products. The draft's lines 160–162 are appropriately cautious; the rewrite should not call the algorithm quadratic.
- **C6 (judgment).** Draft line 56 calls the periodic $\partial_t\Omega$ equation "the Bloch equation [Bloch1958]"; Bloch's equation is time independent, and the periodic form is its Sambe/time-dependent version. The notes cite Mikami, Sadreev, Eckardt–Anisimovas for the Floquet form. Bibliography is the other agent's scope.

## 5. Numerical checks

Run in the `julia` MCP session with `env_path` on this branch, plain `LinearAlgebra` matrices plus one run of the package itself. Exact references come from diagonalizing the truncated Sambe generator $\mathbb G$ ($|m|\le14$) and selecting the $d$ eigenvectors with largest weight in the zero sector: $\Omega_m=V_mV_0^{-1}$, $\mathcal B=V_0\Lambda V_0^{-1}$. The exact van Vleck pair is obtained by solving $\langle\log(\Omega(t)N)\rangle=0$ for $N$ by iteration on a 128-point time grid.

| # | Check | Result |
| --- | --- | --- |
| N1 | Package recurrence vs exact $\mathcal B$ (random generic $3\times3$, harmonics $0,\pm1,\pm2$, seed 1) | $\|\mathcal B-\mathcal B^{[N]}\|$ for $N=0..3$: $\omega=20$: 0.25, 0.043, 0.0086, 0.0016; $\omega=80$: 0.062, 0.0026, 1.3e-4, 5.6e-6; slopes $\omega^{-1..-4}$ |
| N2 | Notes' Riccati (translated) by fixed point vs exact | $\|\mathcal B_{\rm Ric}-\mathcal B\|\le5\times10^{-13}$, $\|X_1-\Omega_1\|\le4\times10^{-15}$ |
| N3 | Trace sum rule | $\operatorname{Tr}\mathcal B^{(1,2,3)}\sim10^{-14}$; $\operatorname{Tr}\mathcal B_{\rm exact}-\operatorname{Tr}\mathcal G_0\sim10^{-13}$; $\operatorname{Tr}\mathcal B_{\rm BM}-\operatorname{Tr}\mathcal G_0\sim10^{-16}$ |
| N4 | Born–Markov | $\mathcal B_{\rm BM}-(\mathcal B^{(0)}+\mathcal B^{(1)}/\omega+\text{2-vertex}/\omega^2)=O(\omega^{-3})$; $\mathcal B_{\rm BM}-\mathcal B$ equals the three-vertex term at $O(\omega^{-2})$ |
| N5 | $\mathcal B^{(2)}$ decomposition | three-harmonic part = van Vleck double sum ($3\times10^{-15}$); two-vertex part $+[\mathcal G_0,N^{(2)}]$ = van Vleck two-vertex sum ($2\times10^{-15}$); gauge term itself has norm 5.3 |
| N6 | Package `BlochFeshbach()` on a qubit Hamiltonian (Pauli, harmonics $0,\pm1,\pm2$, seed 7), `order` 4 | BF = HD components to $10^{-14}$; $H^{(1)}=i\mathcal B^{(1)}$, $H^{(2)}=i(\mathcal B^{(2)}+[\mathcal G_0,N^{(2)}])$ to $10^{-15}$; $\|H_{\rm VV,exact}-H^{[4]}\|\cdot\omega^4\approx84$ at $\omega=120,240,480$; kick harmonics vs exact $i\Lambda_m$: $9\times10^{-4}$, $6\times10^{-5}$, $4\times10^{-6}$ at $\omega=15,30,60$ |
| N7 | Exact pair | $N=\langle e^\Lambda\rangle$ to $10^{-16}$; $\Lambda$ anti-Hermitian to $10^{-13}$; $H_{\rm VV}$ Hermitian |
| N8 | Isometric normalization | see C2 |
| N9 | Brillouin–Wigner roots | $\|\mathcal B_{\rm BW}(\lambda)v-\lambda v\|\sim5\times10^{-14}$ for all eigenpairs of $\mathcal B$ |
| N10 | Rotating-jump RR fixture, $L(t)=\sigma_z+\cos\omega t\,\sigma_y+\sin\omega t\,\sigma_x$ | $L^\dagger L=2$; $R_{\pm1},R_{\pm2}$ only; package formula gives $-i[-\frac{5\gamma^2}{4\omega}\sigma_z,\cdot]$ to $6\times10^{-16}$ (the $+5/4$ sign is off by 5.0) |
| N11 | Self-consistent schemes (item F2) | see Section 6 |
| N12 | Simplex identity (item F1) | see Section 6 |

## 5b. Question 4: the certification record (CP lines 1446–1518)

| Claim | Status | Evidence |
| --- | --- | --- |
| "PR #122 implements" the transported-jump Gram path (1450–1459) | **stale** | PR #122 closed, not merged; code lives only on `origin/feat/105-cp-amplitude-transport`; `git grep` finds no transport code in `src/` on `main` or this branch |
| Inference-clean hot path guarded by JET tests (1460) | **not verifiable** | refers to #122 code; this branch has `test/quality/JET.jl` for its own code only |
| Preserves channel provenance, never calls `positive_completion` (1460) | **not verifiable** | #122 code absent; the provenance on this branch (`microscopic_provenance`, `src/engine.jl:279`) serves positive completion |
| Scope: one-dissipator tranche certified by Sec. oneR (1460) | **not verifiable** (mathematical claim, code absent) | |
| Five-benchmark validation matrix (1463–1471) | **stale as attributed** | similar fixtures exist on this branch as completion tests: driven qubit (`test/cp_completion_validation.jl:153`), full-rank non-diagonal (`:191`, `test/gram_completion.jl:84`), rational nonunitary frame (`:389`), Kerr number-selective loss (`:271`), one- and two-photon loss (`test/cp_completion_analytic_examples.jl:116`); none tests a one-dissipator transport path |
| Nonzero-RR fixture: $H^{\rm VV}_{RR}=-5\gamma^2/(4\omega)\,\sigma_z$, decomposed $1+\frac14$ (1474–1499) | **holds** (mathematics, N10); fixture **absent** from the branch's tests | |
| Rate vs amplitude grading: periodic scalar rates treated as amplitudes are rejected (1502–1508) | **stale / not verifiable** | refers to the #122 path; on this branch `jump(J, γ₀+γ₁cos ωt)` is explicitly supported (`src/liouvillian.jl:404-407`); the odd-onset obstruction exists in completion (`src/completion_types.jl:92-103`, `test/cp_completion_validation.jl:291`) |
| Performance baseline 0.508 ms vs 2.271/8.239/11.958 ms (1511–1515) | **not verifiable** | #122 path absent; no benchmark of it on this branch |
| Coherent solve motivated the fast Bloch–Feshbach backend (1515) | **holds** in the sense that `BlochFeshbach()` now exists (`src/bloch_projection.jl`, `src/expansion_algorithms.jl`, commit 7b96c2c) | |
| Implementation order: retain #122, finish HD and BF canonical algorithms, then open-leg recurrence (1517–1518) | **partly stale** | HD and BF canonical algorithms land in PR #323 (open); #122 is closed, not retained on this stack |
| HD cubic vs BF quadratic bookkeeping (448, outside the record) | **holds for the recurrence only** | engine commutator count $2\sum_n\sum_j(n-j+1)$ is cubic; see C5 for the conversion stage |
| "FE.jl Van Vleck sign and denominator" from the ordered history (560–565) | **holds** | matches the package $\mathcal B^{(1)}$ |

## 6. Question 5: internal correctness of the notes

- **F1 (verified-by-check).** CP eq:charged-first-return-integral (lines 551–558) claims $\int_{0<\tau_1<\dots<\tau_q<2\pi}e^{-i\sum\nu_j\tau_j}=2\pi i^{q-1}\prod_{j<q}s_j^{-1}$ whenever all partial sums $s_j\ne0$ for $j<q$. This fails when two intermediate partial sums coincide: $\nu=(2,-1,1,-2)$, $s=(2,1,2,0)$, exact $0$ against claimed $-i\pi/2$; $\nu=(1,1,-1,-1)$, exact $-3\pi i$ against $-\pi i$; $\nu=(1,2,-2,-1)$, exact $-5\pi i/3$ against $-2\pi i/3$ (exact recursion and independent cumulative quadrature agree). It matched in every tested case with pairwise distinct $s_1,\dots,s_{q-1}$; I did not prove that condition sufficient. The conclusion at line 575 ("This proves that the Fourier partial-sum/fold DAG ... is the correct combinatorial IR") is therefore not established by this identity.
- **F2 (verified-by-check).** NZ line 291 and line 376 state that the self-consistent scheme is the Sylvester pair with $Q\mathbb L_cQ$ neglected, "its only approximation". Dropping $Q\mathbb L_cQ$ leaves the internal line with $\mathcal L_0$: $\mathcal L_0Y_m-Y_m\mathcal L_{\rm eff}-im\omega Y_m=-\mathcal L_m$. The self-energy of eq:nz-self-energy (lines 381–384) instead uses $\operatorname{ad}_{\mathcal L_{\rm eff}}$, dressing the internal line too (as line 197 and line 376's second sentence say). The two generators differ at $O(\omega^{-3})$ ($\|\cdot\|\omega^3\approx8.2$–$8.6$, single-harmonic random example, $\omega=20,40,80$); both miss the exact generator at $O(\omega^{-3})$.
- **F3 (verified-by-check for the premise; likely for the consequence).** NZ footnote at line 193: "For a dissipative $\mathcal L_0$ these [eigenvalues of $\operatorname{ad}_{\mathcal L_0}$] sit in the closed left half-plane". False: $\operatorname{spec}\operatorname{ad}_{\mathcal L_0}=\{\lambda_i-\lambda_j\}$ is symmetric under $i\leftrightarrow j$; a damped qubit gives $\max\operatorname{Re}(\lambda_i-\lambda_j)=+\gamma$. Hence the integrals in eq:nz-markov, eq:nz-ad-resolvent, and eq:nz-graph-Y diverge exponentially on such components, and "the boundary term at infinity removed by ... oscillatory averaging" (line 259) is not a valid step. The resolvent and Sylvester forms are correct as analytic continuations (which the same footnote invokes), and the algebraic invariance equation is unaffected. The phrase "exact Markov closure" should be read as "the algebraic Sylvester problem".
- **F4 (verified-by-check, minor).** NZ line 220: the three-vertex leading term equals the van Vleck double sum "up to the same micromotion gauge". The gauge $[\mathcal L_0,S]$ contains $\mathcal L_0$, so it contributes nothing to the three-harmonic sector; the equality is exact (N5).
- **F5 (verified-by-check for the difference; likely for the literature identification).** CP eq:FeshbachBWVVchain and line 965 identify Bloch ↔ canonical van Vleck with the Shavitt–Redmon/Durand relation. Structurally right (both parametrize the same subspace), but the Floquet van Vleck gauge is a Toeplitz, zero-average single exponent that block-diagonalizes all harmonic sectors, not the P/Q block-off-diagonal canonical normalization; the isometric $S^{-1/2}$ representative differs from the package's by a static unitary at $O(\omega^{-3})$ (C2). The phrase "up to a unitary rotation" (line 200) covers this; the boxed chain does not.
- **F6 (verified).** Mixed Fourier conventions in the CP note (C4) and mixed order indexing (Section 1). No sign error found within any single convention.
- **F7 (uncertain).** NZ convergence theorem (lines 300–308, eq:nz-convergence-condition $\kappa\Lambda\le\frac14$): the contraction bounds check out line by line, but I did not verify the analytic implicit-function step. It concerns the exact-resolvent vertex series, not the package's $\omega^{-n}$ series. Related remark: for bounded finite-dimensional generators the $\omega^{-n}$ series of $\mathcal B$ should also converge above a frequency threshold (analytic perturbation of an isolated Sambe cluster), so the engine docstring's "asymptotic rather than convergent" (`src/engine.jl:219-220`) is the conservative general statement; N1 shows geometric decrease at fixed $\omega$.
- **F8 (likely, minor).** NZ line 324, "their spectra coincide order by order", holds for formal series; finite truncations $\mathcal B^{[N]}$ and $\mathcal G^{[N]}_{\rm eff}$ are not isospectral, only equal to the retained order.
- **F9 (uncertain).** NZ line 278 places the selected branch in the central zone $\operatorname{Im}\lambda\in(-\omega/2,\omega/2]$; at strong drive the branch connected to $Y=0$ need not stay there.
- **F10 (uncertain, not checked).** TCL claims in NZ §5 (quasi-periodic spectral lines at $\omega\pm\epsilon$, finite-time zeros of $\det G(t)$), the phase-blindness and "invariant ring" completeness argument (lines 425–427), and the CP-expansion accuracy statement eq:nz-cp-accuracy.
- Verified without issue: NZ eq:nz-bm-ad and its expansion, eq:nz-three-vertex sign and leading term, eq:nz-trace-sum-rule pairing, eq:nz-gauge-N-expanded, eq:nz-cp-second-order algebra; CP eq:twojump-cumulant-exact, eq:RR-amplitude-pasting, eq:FM1-package and eq:VV1-package, eq:BWresolvent.

## 7. Question 6: what belongs on a public theory page

### (a) Established material for a brief introduction to `BlochFeshbach()` as implemented

- Bloch equation $\partial_t\Omega=\mathcal G\Omega-\Omega\mathcal B$, intermediate normalization $\langle\Omega\rangle=1$, $\mathcal B=\langle\mathcal G\Omega\rangle$ (Bloch; des Cloizeaux for non-Hermiticity; Killingbeck–Jolicard reviews).
- Sambe form $\mathbb G(P+X)=(P+X)\mathcal B$, its P row and Q row; spectrum of $\mathcal B$ = Floquet exponents of the connected branch; non-resonance condition.
- Feshbach/Brillouin–Wigner energy-dependent operator $\mathcal B_{\rm BW}(z)$ and the root relation $\mathcal B_{\rm BW}(\lambda)v=\lambda v$ (Feshbach; Mikami et al.; Eckardt–Anisimovas).
- The $\omega^{-n}$ recurrence and its leading terms, including the explicit $\mathcal B^{(2)}$ of Section 2.3.
- Normalization $\Omega N=e^\Lambda$, $\langle\Lambda\rangle=0$, $N=\langle e^\Lambda\rangle$, $\mathcal K=\Lambda$, $\mathcal G_{\rm eff}=N^{-1}\mathcal BN$, $N^{(2)}$, $\mathcal G^{(2)}_{\rm eff}=\mathcal B^{(2)}+[\mathcal G_0,N^{(2)}]$ (Shavitt–Redmon/Durand for the Bloch vs canonical relation in general; the Floquet zero-average condition is the package's).
- $\operatorname{Tr}\mathcal G^{(n)}_{\rm eff}=\operatorname{Tr}\mathcal B^{(n)}=0$ for $n\ge1$ (Liouville's formula; commutators are traceless).
- For Liouvillians: all maps general, $N^{-1}\mathcal BN$ an ordinary similarity, no CP guarantee; one optional sentence on the NZ reading of the P/Q split.

### (b) CP-native, amplitude-space, and period-map material: keep off the page

NZ §§ 9–15 (frame covariance, CP expansion, physical-average gauge $C^{[N]}$, two-time correlations, dressed channels, Pawula boundary); CP §§ 4–6 (Chruściński–Kossakowski amplitude Feshbach, Floquet lift of amplitudes, Kraus hierarchy, cMPS and HP/SLH representations, open-leg BF/HD recurrences, charged first-return walks, Reimer–Wegewijs half-chains, RR amplitude pasting, finite CPTP period-amplitude triangle); the one-dissipator Gram generator and its static gauge $\mathcal B^{(2)}_R$ (CP §§ 967–1227); Floquet–Davies, HCM/root relation, embeddability. At most one sentence on the page, for example "A construction that expands physical jump amplitudes before assembling the Liouvillian is a separate problem and is not part of this algorithm", with no link until a public page exists (issues #168, #99, #104/#105).

### (c) Open research claims

The NZ reading "the Markov closure is exact" (the invariant-subspace mathematics is established; the NZ integral route is formal, F3); the convergence theorem (F7); ad-dressed resolvent schemes, self-consistent and constrained schemes (F2) and their convergence rates; TCL quasi-periodicity and dressed-projector repair; phase-blindness and invariant-ring completeness; CP-expansion physicality and accuracy $O(\omega^{-(N+1)})+O(\gamma^2/\omega)$; the Lamb-map surjectivity conjecture; the charged first-return identity (F1); the arbitrary-order Floquet–Kraus program; the novelty statements (CP lines 79, 1421).

### Proposed outline, mirroring `high_frequency_expansion.md`

1. **Intro** (mirrors the opening paragraph). One paragraph: the algorithm solves for a wave operator instead of the micromotion, then normalizes to the same van Vleck pair. Names: Bloch, des Cloizeaux, Feshbach, Brillouin–Wigner Floquet expansion.
2. **Wave operator and Bloch equation** (mirrors "Periodic Lie transformation"). Conventions $\dot x=\mathcal Gx$, $e^{-im\omega t}$, $x=\Omega y$, $\dot y=\mathcal By$; the Bloch equation; $\langle\Omega\rangle=1$; $\mathcal B=\langle\mathcal G\Omega\rangle$; for a Hamiltonian $\mathcal B=-iH_B$ with $H_B$ non-Hermitian. Note that literature using $e^{+im\omega t}$ requires $m\to-m$.
3. **Recurrence** (mirrors "Averaging and the homological equation"). $\mathcal R^{(n)}$, $\mathcal B^{(n)}=\langle\mathcal R^{(n)}\rangle$, $\Omega^{(n+1)}_m=\frac im\mathcal R^{(n)}_m$, the same zero-average inverse as the HFE page. State that each order is a linear step with ordinary products, while the full equation is quadratic (fixes C1). Hamiltonian weight $1/m$ in the $H$ language.
4. **Leading terms** (mirrors "Leading terms"). $\mathcal B^{(0)}=\mathcal G_0$, $\mathcal B^{(1)}=i\sum_{m>0}[\mathcal G_{-m},\mathcal G_m]/m$ (equal to van Vleck), and $\mathcal B^{(2)}$ explicitly; trace property.
5. **Sambe-space projection** (mirrors "Sambe-space block diagonalization"). $\mathbb G=-iQ$, $P$ the zero sector; $\mathbb G(P+X)=(P+X)\mathcal B$; P row and Riccati Q row; spectrum of $\mathcal B$; Feshbach/BW $\mathcal B_{\rm BW}(z)$ and the root relation; the recurrence as the expansion of $(\operatorname{ad}_{\mathcal G_0}+im\omega)^{-1}$ in $\operatorname{ad}_{\mathcal G_0}/(m\omega)$, so the package's series is a high-frequency series, not a resolvent (BW) resummation.
6. **Normalization to the van Vleck gauge** (new; the HFE page has no counterpart). $\Omega N=e^\Lambda$, $N=\langle e^\Lambda\rangle$, identification, $N^{(2)}$, $\mathcal G^{(2)}_{\rm eff}$; the three-harmonic part of $\mathcal B^{(2)}$ is already the van Vleck one; Lyndon-basis remark. For a Hamiltonian, the isometric normalization $S^{-1/2}$ gives a Hermitian representative unitarily equivalent to, but in general different from, the van Vleck one (C2).
7. **Liouvillians and complete positivity** (mirrors "Complete positivity after truncation"). General maps, ordinary similarity, no GKSL guarantee for $\mathcal B$ or $\mathcal G_{\rm eff}$; optional NZ sentence (Section 3); link to CP-preserving completion; optional one-sentence pointer from (b).
8. **Comparison with Hori–Deprit.** Keep the current table; state costs as recurrence bookkeeping only (C5).

## 8. Question 7: references the notes rely on for the Bloch/Feshbach and NZ material

As cited in the notes; metadata not verified beyond the notes' own bibliography entries.

- NZ projection and TCL: Nakajima, Prog. Theor. Phys. 20, 948 (1958) [`Nakajima1958`]; Zwanzig, J. Chem. Phys. 33, 1338 (1960) [`Zwanzig1960`]; Breuer and Petruccione, *The Theory of Open Quantum Systems* (OUP) [`Breuer2007`]; Gonzalez-Ballestero, Quantum 8, 1454 (2024) [`GonzalezBallestero2024Tutorial`]; Shibata, Takahashi, Hashitsume, J. Stat. Phys. 17, 171 (1977) [`Shibata1977`]; Chaturvedi and Shibata, Z. Phys. B 35, 297 (1979) [`Chaturvedi1979`]; Ferguson, Zilberberg, Blatter, PRR 3, 023127 (2021); Gu, J. Chem. Phys. 160, 204113 (2024).
- Projector averaging precedent: Buishvili and Menabde, Sov. Phys. JETP 50, 1176 (1979) [`Buishvili1979`, already in `docs/src/refs.bib`]; Saiko, Fiz. Tverd. Tela 35, 20 (1993) [`saikoOne1993`]; Saiko, Theor. Math. Phys. 161, 1567 (2009) [`Saiko2009`].
- Feshbach, Bloch, canonical normalization: Feshbach, Ann. Phys. 5, 357 (1958); 19, 287 (1962); Shavitt and Redmon, J. Chem. Phys. 73, 5711 (1980); Durand and Paidarová, PRA 58, 1867 (1998). The notes do not cite Bloch 1958 or des Cloizeaux 1960 (both already in `refs.bib`).
- Floquet Feshbach / BW / Sambe: Shirley 1965 and Sambe 1973 [`Shirley1965`, `Sambe1973`]; Sadreev, PRE 86, 056211 (2012); Mikami et al., PRB 93, 144307 (2016), erratum PRB 99, 019902 (2019); Eckardt and Anisimovas, NJP 17, 093039 (2015); Wang et al., PRB 110, 245108 (2024) (Floquet Schrieffer–Wolff via Sylvester equations); Keliri and Schirò, arXiv:2606.09727 (2026) (Sambe resolvent, continued fractions); Restrepo et al., PRL 117, 250401 (2016) (transfer tensors).
- Liouvillian van Vleck: Ikeda, Chinzei, Sato, SciPost Phys. Core 4, 033 (2021); Schnell, Denisov, Eckardt, PRB 104, 165414 (2021); Schnell, Eckardt, Denisov, PRB 101, 100301(R) (2020); Venkatraman et al. (static effective Hamiltonians, 2022).
- Amplitude-space (later stacks): Chruściński and Kossakowski, PRL 111, 050402 (2013); Reimer and Wegewijs, SciPost Phys. 7, 012 (2019); Reiter and Sørensen, PRA 85, 032111 (2012); Tokieda et al., PRA 109, 062206 (2024).

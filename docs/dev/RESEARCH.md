# Research record

External findings behind `DESIGN.md`: the SQA API surface, the published oracle, an adversarial
review, algorithms and tooling, and the gauge experiments. Evidence, with sources and confidence.
Decisions drawn from it are in `DECISIONS.md`.


---

## 1. SecondQuantizedAlgebra 0.10.1 — API surface

Source read at `~/.julia/packages/SecondQuantizedAlgebra/DYmwb/`.

### Type hierarchy (NOT the old QSym/QMul/QAdd three-tier)

    QField (abstract)
      QSym (abstract)
        Op            <- the ONLY concrete leaf type
      QAdd            <- the ONLY compound type

- `Op` (`operators/op.jl:33-42`): fields `kind::OpKind`, `name_id::Int32`, `space_index::Int32`,
  `index::Index`, `l1/l2/g/nlev::Int32`. All concrete, `isbits`.
- **`Destroy`, `Create`, `Transition`, `Pauli`, `Spin`, `Position`, `Momentum` are FUNCTIONS**, not
  types. They return differently-tagged `Op`s. Role tests are `is_destroy(o)`, `optype(o)`, never
  `isa`.
- **There is no `QMul` and no `QNumber`.** A monomial is a `QTerm` (`expressions/qterm.jl:14-20`):
  `ops::Vector{Op}`, `ne::Vector{NonEqualPair}`, and a CACHED `hash::UInt`. It is used as a Dict key
  inside `QAdd`.
- `QAdd` (`expressions/qadd.jl:13-19`): `arguments::Dict{QTerm,CNum}` plus `indices::Vector{Index}`.
  So the scalar/operator split lives at the QAdd level: `QTerm` is the noncommutative part, the Dict
  VALUE is the scalar.

### Coefficients: three tiers

`Coeff`/`CNum` (`expressions/cnum.jl:12-16`): `z::ComplexF64` plus
`tail::Union{Native, Poly, Complex{Num}}`.

- `Native` = plain ComplexF64 fast path
- `Poly` = hand-rolled sparse polynomial, kept OFF SymbolicUtils hash-consing for speed
- `Complex{Num}` = full SymbolicUtils fallback

**Promotion is value-dependent** (measured): `(1//2)*a` stays `Native`; `(1//3)*a` falls straight to
`Complex{Num}`, because `Native` is ComplexF64 and 1/3 is not exactly representable. `im*a` and
`2*a` stay `Native`. `(1/3)*a` (Float64) stays `Native`.

Per SQA's own devdocs (`docs/src/devdocs.md:160`): once a coefficient is on the `Complex{Num}` tier,
**coefficient arithmetic dominates runtime, not the operator algebra**.

`to_num(c)::Complex{Num}` (`@public`) is the documented way to read a coefficient at a boundary.

### Normal ordering and simplification

Automatic, not opt-in. `normal_order` docstring (`algebra/algebra.jl:102-131`) verbatim: it is
*"the identity on anything built through public arithmetic: `*, +, -, ^, commutator, Σ, substitute`,
and `adjoint` all canonicalize eagerly"*. So calling `normal_order` in a loop that builds via
`commutator`/`+`/`*` is pure waste.

`simplify` = `normal_order` + `Symbolics.expand` + `SymbolicUtils.simplify` per coefficient. SQA
documents it as a finalizer (`algebra/algebra.jl:179-181`): *"reach for `simplify` as a finalizer…
use `normal_order` for intermediate steps"*.

**Consequence for pruning**: operator-level zeros vanish for free (`_addto_key!`,
`expressions/qterm.jl:170-189`, deletes zero-coefficient entries eagerly, so the dict never stores
one). But coefficient-level cancellation on the `Complex{Num}` tier does NOT happen without
`simplify`. Mathematically-zero buckets therefore survive as structurally nonzero.

### Commutator

`commutator` and `anticommutator` are exported (`SecondQuantizedAlgebra.jl:115`), declared
`algebra/algebra.jl:366`. Always returns `QAdd` (verified). Not rule-based: hand-written
`@enum`-dispatched functions `_can_commute` / `_commute_pair` / `_reduce_pair`
(`operators/operators.jl:367-444`). `commutator(::QAdd,::QAdd)` skips term pairs on disjoint
subspaces as a zero shortcut.

### Hashing and equality — safe

- `Op`: structural `hash`/`isequal` over interned Int32 fields (`operators/op.jl:52-72`). Interning
  (`expressions/intern.jl:14-84`) is idempotent within a session.
- `QTerm`: caches its hash at construction.
- `QAdd`: `isequal`/`hash` delegate to `Dict` isequal/hash, which are **order-independent**
  (`expressions/qadd.jl:128-137`). Verified: `(a+a')*(a+a')` and its expanded form hash equal and
  compare equal.
- `==` returns a plain `Bool` for every type, never a symbolic object.

**So `QAdd` is safe as a Dict key within one session.** Caveat from the source
(`expressions/intern.jl:9-12`, `docs/src/devdocs.md:355`): interned ids are insertion-order within a
session and never serialized. **Never persist a hash-keyed cache to disk.** `Poly` hashes factors on
`objectid`, which is safe only because `@variables` hash-conses.

### Zero — sharp edges

- Zero is the **empty QAdd**. `iszero(a::QAdd) = isempty(a.arguments)` (`expressions/qadd.jl:91`).
- `zero(::Op)` returns the **Julia integer 0**, not an operator zero (`types.jl:30-31`).
- `(a-a) == 0` and `isequal(a-a, 0)` both return **false** silently: there is no
  `==(::QAdd,::Number)` and `QField` is not `<: Number`.
- **`iszero(x)` is the only reliable zero test.**

### Time dependence

Grepped for `fourier`, `harmonic`, `floquet`: **zero hits**. No Fourier/Floquet machinery upstream.
What does exist:

- **`expim(x)`** (`expressions/cnum.jl:39-99`, `@public`): compact unit-phase atom for provably real
  `x`. `expim(x)*expim(y) == expim(x+y)`, `expim(x)*expim(-x) == 1` exactly. Survives `conj`,
  `substitute`, differentiation. **The right primitive for building harmonics.**
- **Bare `exp(im*w*t)` does NOT stay compact.** It becomes `cos(t*w) + im*sin(t*w)`. Attribution
  correction: this is **Base's** `exp(::Complex)` firing before SQA sees anything (leaving vestigial
  `exp(0)` factors), not an SQA behaviour, so it cannot be configured away.
- `exponential_form` / `trigonometric_form` (`algebra/algebra.jl:217-223`, `@public`) convert a
  whole QAdd's coefficients between the two. Verified: `exponential_form(cos(w*t)*a)` gives
  `(0.5exp(im*t*w) + 0.5exp(-im*t*w))*a`.
- **`UnitaryTransform` / `Rotation` / `Displace` / `Squeeze` / `transform` / `conjugate` /
  `gauge_term` / `generators`** (`algebra/unitary*.jl`, exported): exact, closed-form frame changes.
  `Rotation(a, theta, t)` with the `DynamicTime` marker makes `transform` add the `i (dU/dt) U^dag`
  gauge term automatically. This is exactly the spec's KPO preamble (displacement then rotation).
- `average` / `undo_average` / `make_time_dependent` (`average.jl`, exported): QuantumCumulants-style
  moment equations. Downstream consumer, not part of the expansion.
- `normal_to_symmetric` / `symmetric_to_normal` (`algebra/weyl.jl:23-79`): Weyl ordering.

### Numeric backend (exported)

`to_numeric`, `numeric_operator`, `numeric_basis`, `numeric_subbasis`, `numeric_embed`,
`numeric_identity`, `numeric_assemble`, `numeric_assemble_td`, `numeric_materialize`,
`numeric_expect`, `numeric_backend`, `NumericBackend`, `QuantumOpticsBackend`,
`QuantumToolboxBackend`. **Use this for the residual-scaling test; do not hand-roll one.**

### Multi-mode and indexed

`ProductSpace{T<:Tuple}`, `⊗`/`tensor`, auto-detection of `space_index` when the subspace type is
unique. `Index`, `IndexedOperator`, `IndexedVariable`, `Σ`/`∑`, `change_index`, `get_indices`,
`assume_distinct_index`. Bound summation indices live on `QAdd.indices`.

### Performance, measured

`*` on two QAdds distributes eagerly over all `|a|x|b|` term pairs
(`algebra/pipelines.jl:225-250`), each pair fully canonicalized. Measured growth, clean O(n^2) in
summand count: 10 terms 24 us, 20 terms 78 us, 40 terms 280 us, 80 terms 1.17 ms. SQA's devdocs
(`docs/src/devdocs.md:218`) call this intrinsic and by design.

Our own measurement, KPO-scale buckets (~4 terms), 15 harmonics:

    one commutator                     21.9 us
    convolution, Dict-backed  (14x15)   4.32 ms
    convolution, Vector-backed(15x15)   4.36 ms
    210 commutators alone               4.60 ms

**Commutators are ~100% of runtime. Container choice is irrelevant to speed.**

### Type stability, measured

    commutator(::QAdd,::QAdd)             -> QAdd                        concrete
    (::Rational{Int}) * (::QAdd)          -> QAdd                        concrete
    (::Complex{Rational{Int}}) * (::QAdd) -> QAdd                        concrete
    Dict-backed convolution               -> Dict{Int,QAdd}    @inferred OK
    Vector-backed convolution             -> Tuple{Vector{QAdd},Int}     @inferred OK

`QAdd` is concrete; the Union is inside `Coeff.tail`, below anything we return. `(1//3)*(a'*a)`
stays exactly `1//3`.

### Exported vs @public

`@public` but NOT exported, so `@reexport using SecondQuantizedAlgebra` does not forward them:
`QAdd`, `QTerm`, `QSym`, `QField`, `Coeff`/`CNum`, `expim`, `exponential_form`,
`trigonometric_form`, `to_num`, `order_key`, `term_order_key`, `qadd_order_key`, `OpKind`.

`sorted_arguments` is **neither** exported nor `@public`. Use `qadd_order_key`/`term_order_key` for
deterministic enumeration.

**Term iteration IS public, via the Base protocol.** `QAdd` defines `iterate`/`length`/`eltype`
yielding `Pair{QTerm, CNum}`, documented in both the `QAdd` and `QTerm` docstrings ("callers reach
`term.ops` / `term.ne` directly"), and used that way by SQA's own `Base.adjoint(::QAdd)`:

    for (term, coeff) in q

An earlier note in this file and in `DESIGN.md` F4 claimed no public iteration existed and that SQA
needed a new accessor upstream. That was wrong; it generalized from `SymbolicUtils.arguments(q)`
raising `MethodError`. Verified by rerunning the collector prototype with `getfield(q, :arguments)`
replaced by `q`: round-trip still exact.

Operator names must be `Symbol`, never `String`; every constructor guards this.

---

## 2. Oracle — arXiv:2108.02861v2 (Venkatraman, Xiao, Cortiñas, Eickbusch, Devoret)

### Their conventions

| | |
|---|---|
| Fourier | `H(t) = sum_m H_m e^{+i m w t}` (their Eq. 6). **Positive** exponent. |
| static object | `K` (the "Kamiltonian"). **Opposite meaning to the spec's `K`.** |
| generator | `S(t)`, via `e^{L_S}` |
| Lie derivative | `L_S [] = [S,[]]/(i hbar)`. One `1/(i hbar)` per bracket. |
| gauge | `S_avg = 0`, i.e. **van Vleck**, same as the spec. Fig. 3 labels Floquet-Magnus as `S(t_0)=0`. |
| order index | `K^(0) = H_0`, `S^(0) = 0`. Eckardt's `H_F^(n+1)` <-> their `K^(n)`. |
| denominators | `(hbar w)^n` at order n, homogeneous of degree n in the m_i |

### Eq. (8a)/(8b) verbatim

    K^(n)_[k] = { H                                        n = k = 0
                { Sdot^(n+1) + L_{S^(n)} H                  k = 1
                { sum_{m=0}^{n-1} (1/k) L_{S^(n-m)} K^(m)_[k-1]    1 < k <= n+1
                { 0                                        otherwise

    S^(n+1)   = { -int dt osc(H)                                          n = 0
                { -int dt osc( L_{S^(n)} H
                             + sum_{k>1}^{n+1} sum_{m=0}^{n-1} (1/k) L_{S^(n-m)} K^(m)_[k-1] )  n>0

with `osc(f) := f - f_avg`. `K^(n) = sum_k K^(n)_[k]`.

### Translation to the spec's conventions

Substitute `H^V_m = H_{-m}`, then flip every summation index `m_i -> -m_i`, then `hbar -> 1`. All
index arguments are linear homogeneous and all constraints are `L(m) != 0` with `L` linear
homogeneous, so both are invariant under the global flip. Only the denominator changes sign, and it
is homogeneous of degree n. Therefore:

    Heff^(n) = (-1)^n * [ their K^(n), hbar=1, symbols read as our H_m ]

Verified at n=1 against eq:Heff1 and n=2 against eq:Heff2.

**RESOLVED for the generator** (verified at n=1 and n=2). Writing T for the translation itself
(substitute, flip, hbar=1), both rules are clean and neither is graded:

    Heff^(n) =  T( their K^(n) )
    K^(n)    = -T( their S^(n) )

Read against the PRINTED form directly the kick rule looks graded, `(-1)^(n+1)`, because T itself
contributes `(-1)^n`. Always state it as T. Derivation in `R3-recursion.md` §6.

### Term counts (checksums)

    order      0    1    2    3    4    5
    K          1    1    2    8   31  138
    S          -    1    2    6   21   84

Counted independently from the PDF and from their code; they agree.

Coefficient sums (their sign convention): `K^(3)` = -37/24, `K^(4)` = 367/120, `K^(5)` = -4583/720;
`S^(3)` = -8/3, `S^(4)` = 21/4, `S^(5)` = -11.

### Low orders, their convention

    K^(0) = H_0
    K^(1) = [H_m1, H_-m1] / (2 m1 (hbar w))
    K^(2) = [[H_m1,H_0],H_-m1] / (2 m1^2 (hbar w)^2)
          + [[H_m2,H_m1-m2],H_-m1] / (3 m1 m2 (hbar w)^2)

Their summation convention (supplement p.13, verbatim): *"When a Fourier index m_i appears in an
expression it implies the summation over all valid m_i in Z. The (composite) Fourier index of a
Hamiltonian term or that of a commutator … should be non-zero unless it is zero by construction.
The choices of m_i violating this constraint are excluded."*

Order 3 (8 terms) and order 4 (31 terms) transcribed in the agent record; order 5 (138 + 84 terms)
not transcribed. Notable non-generic entries at order 4: term 19 carries `-1/45` with an `(m1-m2)`
denominator, term 21 carries numerator 3 over `40 m1^2 m2 m3`. Both confirmed visually and against
their code, not OCR artifacts.

### Known defects in the paper

1. **`S^(3)`, fifth term, printed denominator `3 m1^3`. Correct is `3 m1^2 m2`.** As printed the
   term vanishes identically, since `sum_{m2} [H_m2, H_-m2]` is odd under `m2 -> -m2`, and nobody
   prints an identically-zero term. Their code emits `3 m1^2 m2`; the order-4 and order-5 analogues
   print `3 m1^3 m2` and `3 m1^4 m2`.
2. **In the printed `S^(5)` block, negative terms appear as a double `+ -` glyph**, not a single
   `-`. It means negative. `K^(5)`, `K^(4)`, `S^(4)` use a single `-`.
3. Page 4 states `e^{L_S} rho = e^{-S/i hbar} rho e^{S/i hbar}`, which is the opposite order to what
   Eq. (2) implies. Anchor on `K = U^dag (H - i hbar d_t) U` with `U = e^{-S/i hbar}`, stated
   consistently on p.2 and supplement p.13.
4. Supplement Eq. (C.6) has a minus where Eq. (4) and (A.1) have a plus. Not load-bearing.

### Worked examples usable as end-to-end tests

**Quantum Kapitza** (their B.9-B.11): `K^(0) = p^2/(2J) - J w_o^2 cos(phi)`, `K^(1) = 0`,
`K^(2) = -(r^2/l^2)(J w^2/8) cos(2 phi)`, `K^(3) = 0`, and `K^(4)` given both normal-ordered (B.10)
and symmetrised (B.11). They note `K^(4)` is NOT obtainable from the classical result by Weyl
quantisation (Groenewold), so the ordering difference is real.

**Driven Duffing** (B.12-B.18): closer to the KPO than Kapitza is. `w_d/w_o = 1.21 ~ 6/5`,
displacement then rotation, `K/hbar = sum_{n=1}^{5} K_n a^dag^n a^n`, with explicit `K_2` and
frequency-shift expansions to `O(1/w^3)`.

**Supplement B III** relaxes "H is order zero in the perturbation parameter" by giving **every node
at row k=0 and column n>=1 a red seed**, i.e. `K^(n)_[0] = H^(n)` for all n rather than `H` at n=0
only. That is the graded-H generalization.

### Reference implementation

`github.com/xiaoxuisaac/vanVleck-recursion` (the URL printed as ref [37]; a commonly-mistyped
variant 404s). Python 3, single 485-line file, sympy only. 3 commits, last 2021. No tests, no
notebooks.

**NO LICENSE FILE. All rights reserved by default.** May be run to generate values; no code and no
vendored output may enter this package.

Their `Term` carries no integer harmonic index at all: only `rotating in {0,1}`, a rational
`factor`, a `freq_denom` integration count, and two children. The `m_i` and the composite
denominators are synthesised at print time by `_latex()`. That is why their output grows `(m1-m2)`
and `(m1-m2)^2` denominators. `integrate()` raises on a DC term.

Verified: the code reproduces Section D exactly (term counts at every order, plus spot checks of
every unusual coefficient).

---

## 3. Adversarial review of the design

### Mathematics: CONFIRMED

The residual and both peeling recursions were re-derived independently from `eq:defining`, then
implemented numerically (dense random 3x3 harmonics, M=2).

- `Heff^(0) = H_0` exact; `Heff^(1)` matches eq:Heff1 to 0.0; `Heff^(2)` matches eq:Heff2 to 8e-15;
  `K^(1)` matches eq:K1 exactly.
- Confirmed by hand: the `-Kdot^(n+1)`-only-at-j=0 split and its index ranges; the `(i/j)` and
  `(i/(j+1))` prefactors; the range `k = 1..n-j+1` (for A from `n-k >= j-1`, for B from
  `n+1-k >= j`, and they coincide); both families vanishing for `j > n`; and that `dressedKdot^(n)_[j>=1]`
  reaches only `dressedKdot^(n-k)_[j-1]` with `k >= 1`, so the order index strictly decreases and the
  recursion is not circular.
- The 1/3: W-coefficients `1/2, 1/2, -1/4, -1/4, -1/6`, giving `-1/4 + 1/2 - 1/6 = 1/12` and
  `1/12 + 1/4 = 1/3`. S-coefficients `1, 1/2, -1/2, -1/2, 0 = 1/2`.
- Reduced convention verified on all three primitives; orders add correctly.

**Residual-scaling validation of orders 3-5**, which the spec has no closed form for. Measuring the
residual of `eq:defining0` against `wd`:

    N        2       3       4       5       6
    slope  -1.037  -1.993  -3.010  -4.005  -5.041

`-(N-1)` rather than `-N` is correct: `eq:defining0` contains `d_t`, which costs one order, and the
omitted `Kdot^(N)` is an order-(N-1) object.

Sharpness demonstrated against three corruption classes: `(i/(j+1)) -> (i/j)`, k-range
`n-j+1 -> n-j`, and a sign flip in `antiderivative` each flatten every slope to -1.

### Errors found in the first design draft

1. **Quasienergy sign.** Written as `Q[m,n] = H_{m-n} + m wd delta`. Under the spec's `e^{-imwd t}`
   convention, `-i d_t e^{-in wd t} = -n wd e^{-in wd t}`, so it is **minus**. The plus belongs to
   the opposite Fourier convention, which would also transpose the off-diagonal to `H_{n-m}`. The
   draft matched neither, and nothing else in the design exercised it.
2. **Missing mean subtraction.** The draft wrote `K^(n+1) = antiderivative(R^(n))` with a comment
   that bucket 0 was "already removed by construction". But `R^(n)_0 = <R^(n)> = Heff^(n)`, nonzero;
   at n=0 it is `H_0`. Correct is `antiderivative(R^(n) - Heff^(n))`.
3. **Two vacuous property tests.** Hermiticity of Heff cannot catch a wrong weight, because `S` and
   `W` are each separately Hermitian, so any real coefficients pass — confirmed to hold even for
   corrupted recursions. `<K^(n)> = 0` is true by construction because `antiderivative` drops
   bucket 0. Neither tests anything.
4. **`normal_order` in the loop is waste**, per SQA's own docstring.
5. **The performance rationale for factoring out wd was wrong** (see the tier measurement in §1).
   The reduced convention is still right, for simplicity.
6. Support bound `(j+1)M` was wrong; `K^(k)` has support `kM`, so the reach is `(n+1)M` independent
   of `j`. The `(n+1)M` conclusion the design used is correct.
7. `Op` is not `<: QAdd`, so `Dict(1 => a)` needs promotion at ingest.
8. The struct had no `wd` and no `t`, yet the API promised "wd reattached" and "K(t) rebuilt".
9. Licensing rule contradicted the test table (which sourced order 5 from the reference impl).
   Also: regenerating oracle values from our own engine is a regression test, not a correctness
   test.
10. `kick(vv)` must exclude `K^(N)`: stage N-1 computes it, but `X^[N] = sum_{k<N}`.

### Review finding that was itself wrong

The review claimed `@inferred`/`@allocations` gates are unachievable for the symbolic core. Measured
and refuted (§1, type stability). `@inferred` applies to the whole engine. Only `@allocations` needs
scoping, since building an SQA expression allocates by construction.

### Review claim not adopted

Last-use freeing of triangle nodes: correctly established that our table cannot be
column-compressed the way a classical Deprit triangle can (every `dressedH^(n)_[j]` reads
`dressedH^(n-k)_[j-1]`
for all `k >= 1`, so nothing can be freed before the final order). But memory is ~720 QAdds at
order 5 and ~5000 at order 10, so it does not matter. Skip unless profiling says otherwise.

### Design gap the review identified

`eq:kickedAgraded` in `sec:comb` is exactly the A-family with the seed replaced by `A_alpha`: same
`(i/j)` peeling, same k range, same convolution. Confirmed out of scope by the user; the seed stays
parameterized regardless.

---

## 4. Algorithms and existing tools

### Deprit vs Hori vs Kamel

Formally equivalent (Henrard, *Equivalence for Lie transforms*, Celest. Mech. 10, 497 (1974);
Mersman). Kamel is Deprit's recursion with machine-convenient bookkeeping. The choice among them is
notation, not complexity.

Only real benchmark found: Nikolaev, arXiv:1612.05207 §8, comparing Gustavson, Hori-Mersman, Deprit,
Henrard, Dragt-Finn and his own in FORM 4.1 to **32nd order** on Hénon-Heiles and Toda-2D.
**Dragt-Finn fastest in every case**; Gustavson slowest. He warns that bracket counts are a poor
proxy because they ignore homological solves, series substitution and memory management.

Koseleff (CeMDA 58, 17 (1994)): Dragt-Finn has the smallest bracket count but larger arguments, and
degrades at very high order from the growing number of exponentiations.

### Improvements over the naive triangle

1. **Nikolaev's summation-free triangle** (arXiv:1612.05207 §6.2):
   `f_k^(n) = f_k^(n+1) + (1/(n+1)) L_{W_{n-k}} f_{n+1}^{(n+1)}`. Same node count, **one bracket per
   node instead of O(n)**. Takes the triangle from O(n^3) to O(n^2) total.
2. **Pymablock's CSE** (arXiv:2404.03728 §3.5): compute `U'^dag X` once and reuse via `+h.c.`,
   factor `A = H'_R U'`, introduce `B = X - H'_R - A`. Effect on Heff: orders 2/3/4 cost 1/3/11
   products versus 1/4/27 for an optimally-implemented reference SW.
3. **Pymablock's coupled Cauchy-product recursion** (§3.3): **O(n) products per added order**.
   Write `U = 1 + W + V`, get `W` from unitarity alone (`W = -U'^dag U'/2`, legal because `U'` has
   no 0th-order term), define `X = [U', H_S]` precisely so `H_0` can be commuted to one side and
   cancelled. Design constraints stated in the paper: one Cauchy product per order; **never multiply
   by H_0** (*"any additional multiplications by H_0 must cancel with additional energy
   denominators… unnecessary work, and it gives longer intermediate expressions"*); exactly one
   Cauchy product by the selected part. `V` comes from one Sylvester solve per order — **free here,
   since our homological operator is diagonal in the harmonic index**.

### The generator: a negative result

**No published algorithm produces a single-exponential generator at better than O(n^2) work per
order.** Pymablock returns `(H_tilde, U, U_adjoint)`; there is no `log` and no generator extraction
anywhere in `pymablock/block_diagonalization.py` (source read). Recovering `K = log(U)` for a
noncommuting graded series costs Theta(n^3) cumulatively by either the associative-powers route or
the dexp^-1 / Bernoulli route, because `d/dlambda log U != U^-1 U'` for noncommuting series.

**One unexplored alternative**: Nikolaev's Kato-resolvent formula (arXiv:1504.05113) gets the
generator directly with no triangle: `F_1 = Q dH/deps`, `F_n = -Q L_V F_{n-1}`, `W = sum_n F_n`,
with `Q` the diagonal small-divisor inverse (our `i/(l w)`). His benchmarks say it wins only for
Hamiltonians with a limited number of perturbation terms and loses at high order for large ones
(Toda-2D, >36000 terms by order 32). ~20 lines; worth one experiment on a driven Kerr oscillator.

### Complexity summary

    naive triangle                O(n^3) total   gives K
    triangle + Nikolaev           O(n^2) total   gives K
    Pymablock                     O(n^2) total   does NOT give K

Asymptotically equal once Nikolaev's variant is in. Only constant factors separate them, including
for the Heff-only case.

### Memory and dependency structure

- A **classical** Deprit triangle depends only on column j-1, so you can sweep column by column and
  keep two columns: O(n) live nodes.
- **Our van Vleck table cannot do this.** `K(n,k)` needs `K(m,k-1)` for all `m < n`, and `S(n)` needs
  the whole row. No node can be freed before the final order.
- Pymablock's answer: instrument the algorithm with a counting object to find entries used exactly
  once and discard them immediately; their code generator emits the deletions.

### Giorgilli / Sansottera — not reusable

Read arXiv:1303.7398 in full. Their algorithm is two lines; the 30+ pages are about **storage**:
packing multi-indices of homogeneous polynomials into a single flat array index via recursively
computed binomial tables with an explicit inverse, plus bit-trie storage for sparse series, plus
(citing Jorba, Experimental Mathematics 8, 155 (1999)) numerical rather than symbolic coefficients.

All of it is for packing **commuting monomials**. Our canonical form is SQA's normal-ordered hash
map. **Do not port.**

Their measured practical ceiling: Hénon-Heiles third integral to **order 58**.

### The series is asymptotic

Giorgilli measures `||Phi_s|| ~ s!`, the theoretically predicted growth, confirmed by root and ratio
criteria. **Optimal truncation order is finite and problem-dependent**: order 32 at moderate energy,
order 9 at high energy, for the same system. Coefficient magnitude growth is intrinsic, not an
implementation defect.

### Existing tools

Full General registry (14 231 packages) extracted and grepped. **No Julia package exists** for
Lie-transform normal forms, Birkhoff normal forms, Deprit, averaging, symbolic Magnus, or symbolic
Floquet expansions.

| Package | What | Noncommuting? | Verdict |
|---|---|---|---|
| NonlinearNormalForm.jl | Dragt-Finn/Lie-map normal forms on GTPSA truncated power series; has `log` of a DA map, `normal()`, `factorise()`. Port of Forest's FPP. | No | Closest Julia prior art, wrong algebra. Its `log`/`exp` map machinery is worth reading. |
| **oameye/VanVleckRecursion.jl** | The user's own. Venkatraman recursion on abstract `Term` objects labelled rotating/static. Generates FORMULAS. Benchmarks at K(5). | structural only | **The other product.** Keep distinct. |
| UnitaryTransformations.jl | Symbolic Schrieffer-Wolff and Lang-Firsov on QuantumAlgebra.jl | Yes | **GPL-3.0.** Hard blocker as a dependency. |
| SecondQuantizedAlgebra.jl | MIT, active | Yes | Correct current dependency. |
| QuantumAlgebra.jl | MIT, 79 stars, active | Yes | Viable alternative backend. |
| Symbolics/SymbolicUtils | Noncommutative symbols **not natively supported** | No | Do not build the operator algebra on it. |
| KOrderPerturbations.jl | DSGE macro models | n/a | Unrelated despite the name. |

**Pymablock** — arXiv:2404.03728, SciPost Phys. Codebases 50 (2025), `quantum-tinkerer/pymablock`,
**BSD-2-Clause**, v2.2.x, active, 95% test coverage. Arbitrary order, multivariate, multi-block,
selective diagonalization. Handles noncommuting symbolic operators via
`pymablock/number_ordered_form.py`, including `LadderOp`/`NumberOperator` with `[a,a^dag]=0`,
`[N,a]=-a`, i.e. Sambe/Floquet ladder operators. **Their tutorial
`docs/source/tutorial/spin_rwa_floquet.md` is our problem**: `H(t) = w_0 sz/2 + g sx cos(Omega t)`
becomes `H_Floquet = w_0 sz/2 + Omega N_a + (g/2) sx (a + a^dag)`, then `block_diagonalize`. Their
Section 5, as of 2024: *"To the best of our knowledge, there are no other packages implementing
arbitrary order quasi-degenerate perturbation theory."*

Their `NumberOrderedForm` keys terms by a tuple of signed powers per mode (negative = creation) with
values that may contain number operators and **non-polynomial functions of them** like `1/(N_b+1)`,
which normal-ordered forms cannot represent at all. Worth understanding before finalizing our
canonical form.

Their `BlockSeries` is a **lazy memoized object**, not a container: `H_tilde[0,0,2]` evaluates if
absent, stores in a dict, returns. Element type is anything with `+` and `*`. NumPy-style slicing
returns a **masked array** so structurally-zero orders and blocks never enter a product. The
optimized algorithm is 14 interlocking series. Because the optimizations would swamp the algorithm's
definition, they wrote a **DSL plus a Python-AST code generator** (`pymablock/algorithms.py`,
`algorithm_parsing.py`). Signal: **the complexity lives in the bookkeeping, not the math.** The
Julia analogue is a macro or `@generated` over a declarative spec.

**SymPT** — arXiv:2412.10240, `qcode-uni-a/sympt`, **MIT**, sympy-based. SWT, full diagonalization,
arbitrary-coupling elimination, operator-level, no Hilbert-space truncation. **Handles time-periodic
Hamiltonians** with the homological equation `[H^0, S^(1)(t)] = -V^(1)(t) + i hbar d_t S^(1)(t)`,
i.e. our operator, and **it does construct S order by order**. Partition-based nested commutators
with caching. Examples to 8th order. Its authors flag their own *"absence of optimization strategies
for reducing the number of commutators"* and point at Pymablock.

**Celestial mechanics**: `celmech` (Hadden, GPL-3.0) exposes only `FirstOrderGeneratingFunction`,
first order in the mass ratio — not an arbitrary-order tool. TRIP (Gastineau & Laskar, IMCCE) is a
mature closed-source Poisson-series CAS, not a library, commuting variables.
**"XSeries" could not be verified to exist**; what exists under Giorgilli's name is a 1979 FORTRAN
CPC program and *Norma*, a REDUCE package.

**BCH**: `bch` by Casas & Murua (arXiv:2102.06570), C, free Lie algebra in Lyndon/Hall/right-normed
bases, exact rationals, hash tables. Degree 20 -> 111 013 coefficients in <0.5 s, 11 MB; degree 30 ->
74 248 451 coefficients, 55 hours, 5.5 GB. The right calibration if a BCH recombination is ever
needed. No Julia equivalent.

### Practical order ceilings

| setting | order |
|---|---|
| commuting, flat packed arrays (Giorgilli) | 58 (32 useful before divergence) |
| commuting, FORM (Nikolaev) | 32 |
| noncommuting symbolic (Pymablock, SymPT) | 8 |
| abstract terms (VanVleckRecursion.jl) | 5 |

**The gap is term proliferation in the noncommutative algebra, not the recursion.**

### Arithmetic traps

- Rational denominators accumulate `k!`-like factors from the `1/k` at every bracket depth times
  Fourier index products. Rough estimate at order 8 with M=8: index products ~72^8 ~ 7e14 times
  Deprit denominators ~8! ~ 4e4, i.e. ~3e19, past `Rational{Int64}`'s 9.2e18. Julia throws
  `OverflowError` rather than wrapping. **Plan `Rational{BigInt}` (or Int128 with a documented
  ceiling) from the start**; retrofitting a coefficient type through a symbolic expression tree
  later is miserable.
- Intermediate expression swell dominates final size; the normal-form literature reports
  intermediates occupying gigabytes while final coefficients print compactly. **Canonicalize at
  every node, do not defer.**
- Jorba's counter-position: drop exact rationals, carry numerical coefficients, accept numbers
  instead of formulas. Worth offering as a mode past order ~10.


---

## 5. Gauge experiment — van Vleck vs des Cloizeaux

Run to discharge the two-engine prerequisite. It produced a negative result.

Setup: random Hermitian harmonics, d=3, |m| <= 2, `H_{-m} = H_m^dag`. Compared the Deprit-triangle
van Vleck `Heff^[N]` against an exact canonical block diagonalization of the truncated Sambe
operator `Q[m,n] = H_{m-n} - m*wd*delta`, via the des Cloizeaux direct rotation
`U = X (X^dag X)^{-1/2}` with `X = Pt P + (1-Pt)(1-P)`, `Pt` the spectral projector onto the band
connected to block 0. Operators compared, not spectra, since spectra are gauge invariant.

    N     slope     expect
    1    -1.000      -1      agree
    2    -2.005      -2      agree
    3    -2.986      -3      agree
    4    -3.020      -4      DISAGREE
    5    -3.017      -5      DISAGREE

Not a truncation artifact: identical at Mcut = 14, 24 and 40, and N=4, N=5 give the same residual
(4.27e-08 vs 4.34e-08), i.e. adding order-4 terms does not reduce it.

Diagnostics on the des Cloizeaux transformation:
- `||Heff(block 0) - Heff(block 1)|| ~ 1e-13` — translation covariant, as Floquet symmetry demands.
- Generator `G = log U` is **not Toeplitz**: deviations among blocks with equal `m-n` are 4.0e-2 at
  |l|=1 and 3.1e-2 at |l|=2, against `||G|| = 5.8e-2`. It is the transformation for one target
  block, not the global periodic K(t).

**Conclusion.** van Vleck and des Cloizeaux are distinct canonical schemes (this is standard: van
Vleck, des Cloizeaux and Bloch all differ; the first two are both Hermitian and agree only at low
order). They coincide through order 2 and diverge at order 3. So "canonical block diagonalization"
does not uniquely determine Heff, and the assumption that Pymablock's gauge equals van Vleck's
cannot be discharged by general reasoning. It must be checked against their specific construction,
at order 3 or higher, because order 2 provably cannot discriminate.

Two reusable tests came out of this: the order-scaling operator comparison, and the
generator-Toeplitz check (`G[m,n]` depending only on `m-n` is necessary for a van Vleck generator).

### des Cloizeaux is symbolic, which raises the risk rather than lowering it

des Cloizeaux is NOT inherently numerical. It is the Bloch wave operator (standard order-by-order
recursion via the generalized Bloch equation) followed by a symmetric orthonormalization,

    H_dC = (P Om^dag Om P)^{-1/2} . P Om^dag H Om P . (P Om^dag Om P)^{-1/2}

with `(1+X)^{-1/2} = 1 - X/2 + 3X^2/8 - ...` terminating at each order because X starts at second
order. The eigen-decomposition route above was used only because it is the cheapest EXACT reference.

Consequence: des Cloizeaux is a legitimate rival engine, not merely a numerical yardstick. And
Pymablock's `W = -U'^dag U'/2` from unitarity is structurally a symmetric orthonormalization, i.e.
it carries the des Cloizeaux signature. This makes the gauge risk more likely, not less.

### FAILED experiment, recorded so it is not repeated

Hypothesis: the discriminator is not des-Cloizeaux-vs-van-Vleck but whether the block partition is
translation covariant, since the single-block projector P breaks Sambe translation symmetry by
construction. Test attempted: take `G = log U_dC`, average blocks along each diagonal to force
Toeplitz structure, re-antihermitize, re-exponentiate, extract Heff.

Result: slope -1.000 at every N from 2 to 5, i.e. it fails to reproduce even order 1.

Diagnosis: the naive diagonal-averaging projection destroys the block-diagonalization property, so
the projected operator is not a valid transformation and its Heff is meaningless. **The experiment
tested the hack, not the hypothesis.** The hypothesis remains untested and plausible.

Two properties must be separated when anyone retries this:
  (a) block-off-diagonal (no m=n blocks)      <- the Schrieffer-Wolff / van Vleck condition
  (b) Toeplitz (depends only on m-n)          <- the Floquet periodicity condition
van Vleck requires BOTH. The experiment above measured only (b). A correct retry must construct a
transformation that satisfies both by construction and still block-diagonalizes exactly, not
post-process one that does not.


---

## 6. Prior art and algorithm choice

Full General registry (14 231 packages) grepped: there is **no** Julia package for Lie-transform
normal forms, Birkhoff normal forms, Deprit, symbolic Magnus, or symbolic Floquet expansions.
Nearest neighbours: `UnitaryTransformations.jl` (symbolic Schrieffer-Wolff, **GPL-3.0**, hard
blocker as a dependency), `NonlinearNormalForm.jl` (Dragt-Finn on commuting phase-space variables,
wrong algebra). Symbolics.jl does not support noncommutative symbols natively. Staying on
SecondQuantizedAlgebra (MIT, active) is correct.

### Position relative to VanVleckRecursion.jl

`oameye/VanVleckRecursion.jl` computes the recursion at the level of *abstract term structures* and
emits formulas. It benchmarks at order 5, because the formula grows combinatorially with nothing to
collapse it.

**FloquetExpansions.jl is the other product**: evaluate the recursion on concrete SQA operators, so
the algebra collapses terms at every node. That is the same nested-sums-vs-sequential-convolutions
argument as R3 notes §6. Keep the two separate and deliberate; do not merge them.

### Is Deprit still the right algorithm?

Deprit / Hori / Kamel are formally equivalent (Henrard 1974, Mersman); the choice among them is
notation, not complexity. Two things do beat the naive triangle:

**Pymablock** (arXiv:2404.03728, SciPost Phys. Codebases 50, BSD-2-Clause) achieves **O(n) operator
products per added order** where the triangle costs O(n^2). Mechanism: write U = 1 + W + V, obtain
W from unitarity alone (W = -U'^dag U'/2), and define the auxiliary X = [U', H_S] so that H_0 is
never multiplied. Their Floquet tutorial is literally our problem in Sambe space, with the Sylvester
solve being our free (i/(l*w)) division.

**But it does not produce the generator.** `block_diagonalize` returns `(H_tilde, U, U_adjoint)`;
there is no log and no generator extraction in the source. Recovering K = log(U) for a noncommuting
graded series costs Theta(n^3) cumulatively by either the power-series or the dexp^-1 route, i.e.
exactly what the triangle costs. Since a materialized K(t) is a hard requirement here, **the triangle
is not beaten for our requirement**. Their CSE work is still worth reading: orders 2/3/4 need
1/3/11 products for them versus 1/4/27 for an optimally-implemented reference.

**Nikolaev (arXiv:1612.05207) — PAPER READ, and it is not what the summary said.**

Earlier drafts of this document called §6.2 a "summation-free triangle" that would drop our
triangle from O(n^3) to O(n^2). That was a misreading. §6.2 computes the TRANSFORM given W; it does
not produce W. Nikolaev is explicit: *"we compute the generator independently and can use faster
algorithms."* His §6.2 recursion

    f_k^(n) = f_k^(n+1) + (1/(n+1)) L_{W_{n-k}} f_{n+1}^{(n+1)},    k = 0..n
    seed f_n^(N) = sum_{k=n}^{N} eps^k H_{k-n},   answer H~ = f_0^(0)

is one bracket per node over O(N^2) nodes, but every `W_{n-k}` on the right must already be known.
It cannot replace a triangle that determines K and Heff interleaved, which is what ours does.

**The real prize is §6.1: an explicit, NON-recursive generator.**

    F_1 = Q  dH/deps
    F_n = -Q  L_V F_{n-1},      n = 1..N            V = sum_{k>=1} eps^k H_k
    W   = sum_{n=1}^{N} F_n[z^n]

Q is the reduced resolvent: it annihilates the kernel of the homological operator and divides by
the small divisor off it. **In our Fourier setting Q IS `antiderivative`** — diagonal in the
harmonic index, `i/(l*wd)` on l != 0. His `z` bookkeeping exists only to avoid negative powers
during CAS manipulation; implementing Q directly as the diagonal division should remove it.

This is a LINEAR CHAIN, not a triangle: N steps, each one Cauchy product plus a diagonal division.

    ours (Deprit triangle)   O(n^2) nodes x O(n) brackets  =  O(N^3) convolutions
    Nikolaev 6.1 + 6.2       O(N^2) + O(N^2)               =  O(N^2) convolutions

A genuine one-power improvement, and it produces the single-exponential generator directly, which
is the thing Pymablock cannot do.

**Gauge: UNRESOLVED. An earlier claim here that it is "NOT van Vleck, do not adopt" is WITHDRAWN.**

Mapping, VERIFIED (`docs/dev/experiments/nikolaev_experiment.jl`). Set eps = 1/wd and divide the Sambe operator by
wd: `Q/wd = H_0 + eps*H_S` with `H_0 := -i d_t / wd` acting as multiplication by `-l` on harmonic
l. Then `L_F G = -i[F,G]`, `L_{H0} X|_l = i*l*X_l`, so `P` = bucket 0 and `S|_l = 1/(i*l)`, which is
`antiderivative` up to a sign. `H_1 = H_S`, all higher `H_k = 0`.

    Heff^(0)   ||triangle - Nikolaev|| = 0.000e+00        VERIFIED
    Heff^(1)   ||triangle - Nikolaev|| = 0.000e+00        VERIFIED
    Heff^(2)   ||triangle - Nikolaev|| = 2.236e-15        VERIFIED
    Heff^(3)   ||triangle - Nikolaev|| = 1.306e+00        UNTRUSTED, see below

Generator, `K^(n+1) = -W^(n)/(n+1)` (the 1/(n+1) is Deprit's dU/deps convention, since the
single-exponential generator is the eps-integral of his W): residual 0 at n=0, 2.6e-16 at n=1.

**Why order 3 is untrusted.** A spectral test settles it: conjugation preserves spectrum, so if the
two differed only by gauge they would be isospectral. Measured at N=4: operator difference 1.31,
**spectrum difference 0.62** against ||H|| = 11.6. So the order-3 discrepancy is NOT a gauge.

But `Heff`'s spectrum IS the physical quasienergy spectrum, which any correct normal form must
reproduce in any gauge. A method published and benchmarked to 32nd order does not get quasienergies
wrong. The error is therefore far more likely in the transcription of his eps^4 `H~` line, which is
eight terms dense with `S^2`/`S^3` superscripts, than in his method.

**Also withdrawn: the inference from `<W^(2)> != 0`.** His `W` is Deprit's dU/deps generator, not
`log U`. The van Vleck condition is on `log U`. `<W> != 0` therefore proves nothing about the gauge.

**Gauge conversion is CHEAP, contrary to an earlier claim here.** If two schemes differ only by
gauge, `D = U_1^dag U_2` is block-diagonal AND Toeplitz (Nikolaev's construction uses only H_S, S
and P, all Toeplitz-preserving). Block-diagonal + Toeplitz forces D to be a STATIC operator. So:

    Heff_vanVleck = e^{-d} Heff_other e^{d}        a static conjugation, cheap
    d solved order by order from the <log U> = 0 condition, with BCH in which one argument is
    always bucket-0, hence far cheaper than a general BCH

Nothing here costs Theta(n^3). So "express an O(N^2) scheme in the van Vleck gauge" is a reasonable
plan, not a dead end.

**PARKED. Status after term-by-term verification** (`docs/dev/experiments/nikolaev_experiment.jl`):

    VERIFIED    eps=1/wd, H_0 := -i d_t/wd mapping
    VERIFIED    Heff^(0), Heff^(1), Heff^(2)  exact to machine precision
    VERIFIED    K^(n+1) = -W^(n)/(n+1)  at n=0, n=1
    VERIFIED    transcription of all 8 H_1-only terms of the eps^4 line, checked term by term
                against the rendered page 12 (composition order, S powers, coefficients)
    VERIFIED    his identity  P L_F P F == 0    holds exactly in our mapping
    FAILS       his identity  P L_F S F == 0    gives 6.565, and 6.565 = 2 x ||Heff^(1)||
    FAILS       Heff^(3): operator differs by 20%, and SPECTRA differ by 0.62 vs ||H||=11.6

The spectral failure is the decisive one: `Heff`'s spectrum is the physical quasienergy spectrum
and is gauge invariant, so a correct scheme must reproduce it in any gauge. A method published and
benchmarked to 32nd order does not get quasienergies wrong.

The picture is therefore internally inconsistent: the mapping reproduces his results exactly at
three orders AND one of his two identities, while violating the other identity and order 3. Read
literally, `P L_F S F == 0` would force `Heff^(1) = 0`, contradicting eq:Heff1. So a notational
subtlety in his formalism is uncaptured — most plausibly the distinction between the PERTURBED
`S_H` (which appears in his defining `W = S_H dH/deps`) and the UNPERTURBED `S_{H0}` that the
displayed expansions suppress the index on.

**Decision: park.** Resolving it means reverse-engineering a terse classical-mechanics formalism
for a speedup on an optimization path that blocks nothing. The triangle is verified and v1 is
unblocked.

**To resume**, start here and do not redo the above:
1. Settle the `S_H` vs `S_{H0}` question in §4/§5 of the paper. That single distinction most likely
   explains both the failed identity and order 3.
2. Gate every attempt on ISOSPECTRALITY, not on operator equality. A correct transcription must be
   isospectral with the triangle at every order; only then is it meaningful to ask whether the
   operators also match (same gauge) or merely the spectra (gauge difference, convertible).
3. If it is a gauge difference, the conversion is cheap: D = U_1^dag U_2 is block-diagonal AND
   Toeplitz, hence a STATIC operator, so `Heff_vanVleck = e^{-d} Heff_other e^{d}` and the BCH
   solving for d always has a bucket-0 argument. Nothing costs Theta(n^3).

Until then the triangle stays the plan of record: it is verified, and it is O(N^3) rather than
wrong.

Not reusable: the Giorgilli/Sansottera school (arXiv:1303.7398). Their entire optimization budget
goes to packing *commuting* monomials into flat arrays via indexing functions and bit-tries. Our
canonical form is SQA's normal-ordered hash map. Do not port their machinery.

### Practical order ceiling, calibrated

| setting | order reached |
|---|---|
| commuting, flat packed arrays (Giorgilli) | 58 |
| noncommuting symbolic (Pymablock, SymPT) | 8 |
| abstract terms (VanVleckRecursion.jl) | 5 |

The gap is term proliferation in the noncommutative algebra, not the recursion. "Arbitrary order"
in practice means roughly 6-10 for a rich drive. Size expectations accordingly.

Also worth reading, MIT: **SymPT** (arXiv:2412.10240) constructs S order by order for time-periodic
Hamiltonians using exactly our homological equation. Its authors flag their own lack of
commutator-reduction optimization and point at Pymablock.

## 7. SecondQuantizedAlgebra facilities we should use rather than rebuild

R1's survey turned up several things the earlier drafts of this document ignored.

**Numeric backend, for the residual-scaling test.** SQA exports `to_numeric`, `numeric_operator`,
`numeric_basis`, `numeric_embed`, `numeric_assemble`, `numeric_assemble_td`, `numeric_expect`, and
backends `QuantumOpticsBackend` / `QuantumToolboxBackend`. §7's primary gate calls for "a
finite-dimensional numeric backend"; do not hand-roll one, use this.

**Frame transformations, for the input pipeline.** `UnitaryTransform`, `Rotation`, `Displace`,
`Squeeze`, `transform`, `conjugate`, `gauge_term`, `generators`. Critically `Rotation(a, theta, t)`
with the `DynamicTime` marker: `transform` then adds the `i (dU/dt) U^dag` gauge term
automatically. The spec's own KPO example requires exactly a displacement followed by a rotation
before the expansion starts (eq:kpo-zeta, eq:kpo-disp, eq:kpo-hs). This is the natural input
pipeline and it already exists, exactly and non-perturbatively. Document it as the recommended
preamble; do not reimplement it.

**Weyl ordering, for one oracle.** `normal_to_symmetric` / `symmetric_to_normal`. Venkatraman B.11
gives the quantum Kapitza `K^(4)` in SYMMETRISED form (B.10 is the normal-ordered equivalent). Use
the converter rather than transcribing both. They also note K^(4) is not obtainable from the
classical result by Weyl quantisation (Groenewold), so this is a real ordering difference, not a
cosmetic one.

**Multi-mode and indexed systems: take a stance.** `ProductSpace`, `tensor`/`⊗`, `Index`,
`IndexedOperator`, `Σ`/`∑`, `change_index`, `assume_distinct_index`. Nothing in the recursion cares
how many modes there are, so multi-mode should work unchanged; the harmonic buckets are still
QAdds. Indexed sums are less obvious: `Σ` bodies interact with normal ordering and with
`QAdd.indices`, and the collector must not lose bound indices when it splits off phase factors.
Decide whether v1 supports indexed sums or explicitly rejects them with a clear error. An
untested "probably works" is the worst option.

**Downstream interop.** `average` / `undo_average` / `make_time_dependent` turn an operator
expression into moment equations. That is the natural consumer of `effective_hamiltonian(vv)` and
worth one line in the docs, not code.

Also from R1: operator names must be `Symbol`, never `String` (every constructor guards this);
`sorted_arguments` is neither exported nor `@public`, so use `qadd_order_key`/`term_order_key` for
deterministic enumeration.


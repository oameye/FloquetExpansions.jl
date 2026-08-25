# FloquetExpansions.jl — van Vleck expansion: design

What to build. Decisions and their rationale are in `DECISIONS.md`; the derivation and its
verification in `DERIVATION.md`; external findings (SQA API, the published oracle, prior art) in
`RESEARCH.md`.

Spec: `tmp/floquet_bath_tracing.tex` (not tracked) sec:vv. Scope: van Vleck gauge only, fully symbolic, on
SecondQuantizedAlgebra 0.10.

## 1. Representation: factor out wd by order

Every object at order n carries exactly wd^-n. Strip it and track the power in the order index.
Under Xred^(n) := X^(n) * wd^n:

    Kdot^(k) wd^(k-1)      = -i l * (K^(k) wd^k)
    antiderivative(X) wd^(n+1)   = (i/l) * (X wd^n)
    [K^(k), dressedH^(n-k)] wd^n     = [K^(k) wd^k, dressedH^(n-k) wd^(n-k)]     orders add

So wd never appears in the *bookkeeping*. The rationale is simplicity: no symbolic wd compounding
through every denominator, and one less thing to get wrong at the order boundary.

Note this is not a speed argument. SQA's `Native` coefficient tier is `ComplexF64`, so promotion is
value-dependent: `1//2` stays `Native`, but `1//3` falls straight to the `Complex{Num}` tier. The
Deprit weights are 1/j and 1/(j+1), so 1/3 appears at `dressedKdot^(n)_[2]` and essentially every
node from order 2 up is on the slow tier regardless. See decision D5.

Caveat: this concerns the recursion's own bookkeeping. A seed harmonic may itself carry wd
(the KPO's xi = -2iF wd/(w0^2-4wd^2)). That is fine as long as it stays an opaque symbol.

## 2. Core type

```julia
struct PeriodicOperator
    components::Dict{Int, QAdd}   # l => X_l
end
```

`Dict` for uniform handling of gapped and dense drives. (It is NOT a sparsity win on the flagship
case: the KPO has H_0..H_7 all nonzero, so support is a contiguous interval at every order. Int
hashing is negligible beside the commutators, so Dict is still the right call, just not for the
originally stated reason.)

Order is not a field; it lives in the memo key.

Ingest promotes: a bare `Op` is not `<: QAdd`, so `Dict(1 => a)` needs `normal_order(a)` or `1*a`.
The constructor does this.

Zero: use SQA's `iszero`, never `== 0`. SQA defines no `==(::QAdd, ::Number)`, so `x == 0` is
silently `false`. `time_average` returns `zero(SQA.QAdd)` when bucket 0 is absent, not `0`
(`zero(::Op)` returns the integer 0).

### The four operations

    commutator(K, X)_l       = sum_p [K_p, X_{l-p}]
    derivative(X)_l          = -i*l * X_l
    time_average(X)          = X_0                     :: QAdd
    antiderivative(X, gauge) = (i/l) * X_l  for l != 0, plus the gauge's constant at l = 0

`time_average`, not `average`: SQA EXPORTS `average` for a different operation, the expectation
value (`average(a'*a)` is `⟨a' * a⟩`), and we `@reexport` SQA, so that name is already in every
user's scope meaning something else. `commutator` by contrast IS the same operation on a new type,
so extend `SQA.commutator` rather than shadow it.

The other names are the spec's. Its sec:vv paragraph on eq:recursion-split calls the K equation the
**homological equation**, `<.>` the **period average**, and names the inverse in so many words:
"its *zero-average antiderivative* is `sum_{l!=0} (i/(l wd)) X_l e^{-i l wd t}`". We keep
`antiderivative` and promote the "zero-average" qualifier to the gauge argument, because that
qualifier is precisely what Floquet-Magnus changes.

The three non-bracket operations are the homological triple, and that is why they always appear
together as `antiderivative(X - time_average(X), gauge)`:

    derivative       the homological operator, the thing the recursion must invert
    time_average     projector onto its kernel (constants are the only periodic functions with
                     zero derivative), which is why <R> is unreachable by any kick and survives
                     as Heff
    antiderivative   inverse on the complement of that kernel, unique only up to a constant

**That leftover constant IS the gauge, so it is an argument, not a default:**

    <K> = 0      -> constant is 0, bucket 0 simply dropped     van Vleck   (v1)
    K(t_0) = 0   -> constant fixed by evaluating at t_0        Floquet-Magnus

Same engine, same triangle, one line different. Dispatching on the gauge type from the start makes
Floquet-Magnus an added method rather than a refactor. See `DECISIONS.md` D7.

`antiderivative` requires a mean-zero input and asserts it. Callers pass `X - time_average(X)`. This is
not a convention, it is the operator's domain: it inverts d/dt, whose image is exactly the mean-zero
periodic operators. A DC term reaching it is a bug by definition. The assert is free and catches
sign errors at the order they are introduced.

Rejected: `antiderivative0` (the `0` was the van Vleck gauge encoded as a digit, and becomes a lie
under Floquet-Magnus); `solve_homological` (right theory word, but does not say which direction the
equation is solved); `integrate` and `primitive` (the first collides with `time_average`, which
is itself a time integral; the second with "the four primitives" right above).

## 3. Collector

    harmonics(H::QAdd, w, t) -> PeriodicOperator     parse symbolic time dependence
    (X::PeriodicOperator)(w, t) -> QAdd              rebuild it; exact inverse
    PeriodicOperator(Dict(m => H_m, ...))            direct

IMPLEMENTED. Two things the implementation added beyond the prototype:

- **Constant phase offsets are supported.** `cos(w*t + phi)` is an ordinary drive, so the parser
  splits the exponent as `c*w*t + offset` (by substituting `w,t => 1` and `w,t => 0` and
  differencing), takes `m = -c`, and returns the offset as an `expim(phi)` factor on the
  COEFFICIENT. It also verifies `arg == c*w*t + offset` symbolically, which is what rejects
  `w*t^2` and anything else not linear in `w*t`.
- **The offset phase must meet the QAdd, not the coefficient.** `expim` returns an SQA `Coeff`,
  and neither `BasicSymbolic * Coeff` nor `BasicSymbolic * Complex{Num}` has a method (the latter
  hits a Symbolics bug in `unwrap(::Complex{Num})`). So `_phase_of` returns the offset separately
  and `harmonics` applies it as `contribution * expim(offset)`, where SQA's own `QAdd * Coeff`
  handles it.

Trap for anyone touching this file: `iszero` on a `BasicSymbolic` builds the symbolic equation
`0 == 0` instead of returning a `Bool`, so it is unusable in a condition and silently poisons a
ternary. Compare structurally against the literal instead.

Parsing:
- `exp(im*w*t)` is expanded to cos+i*sin by **Base**'s `exp(::Complex)`, before SQA ever sees it
  (leaving vestigial `exp(0)` factors that pin the coefficient to the slow tier). Not an SQA
  behaviour and not configurable. Normalize input through `exponential_form` first; verified that
  `exponential_form(cos(w*t)*a)` recovers `(0.5exp(im*t*w) + 0.5exp(-im*t*w))*a`.
- `expim(x)` is the compact phase atom, and the right primitive for *building* phases on output.
- `H_{-m} = H_m^dag` (eq:fourierH) is CHECKED, not enforced here: `harmonics` is a pure parser and
  the coupling operator will run through it too. `ishermitian(::PeriodicOperator)` exists for it,
  and `floquet_expansion` is the ingest point that must assert it, so a malformed drive still
  fails before it becomes a non-Hermitian Heff at order 3.

### Index extraction — PROTOTYPED AND WORKING

Working prototype at `docs/dev/experiments/collector_prototype.jl`, round-trip verified exact
(`iszero(simplify(rebuilt - H))`).

**F1. Coefficients merge across harmonics. This breaks the naive one-harmonic-per-term design.**
`QAdd` is a map monomial -> coefficient, so two input terms with the SAME monomial but different
phases merge into ONE dict entry whose coefficient is a SUM over harmonics. Observed:

    input :  2*expim(-w*t)*a  +  5*expim(2*w*t)*a
    stored:  (2exp(-im*t*w) + 5exp(im*2t*w)) * a        <- one entry, two harmonics

So the collector MUST split each coefficient into additive parts and classify each part
separately. A per-term harmonic lookup is wrong.

**F2. The operation is `expim`, and its argument is the REAL exponent.** It DISPLAYS as
`exp(im*2t*w)` but `SymbolicUtils.operation` returns `expim`, and `arguments(x)[1]` is `2t*w`, not
`im*2t*w`. So recover the exponent by substituting `w => 1, t => 1`, which is far simpler than the
`w => -im` trick. Then, for our `e^{-i m w t}` convention:

    expim(c*w*t) == exp(+i c w t) == exp(-i(-c) w t)     =>   m = -c

**F3. The storage trap does not surface here.** R1 warned that a negative phase is stored as the
positive atom at exponent -1. Through the `to_num` view it comes back as `expim(-t*w)` directly,
no power. The `numerator(exponent)` handling is unnecessary for this route.

**F4. RETRACTED. Public term iteration exists; no upstream change is needed.** An earlier note
here claimed the collector had to reach `getfield(q, :arguments)` and that SQA needed a new public
accessor. Wrong on both counts. `QAdd` implements the Base iteration protocol:

    Base.iterate(q::QAdd)     -> Pair{QTerm, CNum}
    Base.length(q::QAdd)
    Base.eltype(::Type{QAdd}) = Pair{QTerm, CNum}

and the `QAdd` docstring specifies it ("Iterating over a `QAdd` yields `Pair{QTerm, CNum}`
entries"), as does `QTerm`'s ("callers reach `term.ops` / `term.ne` directly"). SQA's own
`Base.adjoint(::QAdd)` is written `for (t, c) in q`. So the collector loop is

    for (term, coeff) in q          # public, documented, no getfield

Verified: the prototype round-trips exactly with the `getfield` replaced by `q`.
What misled the survey was that `SymbolicUtils.arguments(q)` raises `MethodError` and
`sorted_arguments` is neither exported nor `@public` — both true, and both irrelevant, because
iteration is the documented route and `qadd_order_key` covers deterministic ORDER when it matters.
Everything else the collector touches is public too: `to_num`, `expim`, `Op` arithmetic,
`one(QAdd)`.

Sketch (full version in the prototype):

    addparts(x)      -> arguments if op is (+), else [x]
    phase_of(part)   -> (m, rest) by finding the single expim factor among a Mul's arguments
    harmonics(q,w,t) -> for each (monomial, coeff): to_num, split real and imag, split additive
                        parts, classify, accumulate rest*monomial into bucket m

### Quasienergy matrix

    QuasienergyOperator(hs, nmax) -> Q,   Q[m,n] = H_{m-n} - m*wd*delta_mn

MINUS on the diagonal. With H_S(t) = sum_m H_m e^{-i m wd t} and u(t) = sum_n u_n e^{-i n wd t},
the operator Q = H_S - i*d_t gives -i*d_t e^{-inwd t} = -n wd e^{-inwd t}. The +m wd diagonal
belongs to the opposite Fourier convention, which would also transpose the off-diagonal to H_{n-m}.

Inspection and interop view only; the recursion reads the harmonic dict directly. Because nothing
else exercises it, it needs its own test (§8).

## 4. Engine: the Deprit triangle

```julia
struct FloquetExpansion{M}
    H           :: PeriodicOperator
    K           :: Vector{PeriodicOperator}  # K[k]      = K^(k),      k = 1..N
    Kdot        :: Vector{PeriodicOperator}  # Kdot[k]   = d/dt K^(k)
    dressedH    :: Vector{PeriodicOperator}  # (i ad_K)^j/j!     on H,    flat over (n,j)
    dressedKdot :: Vector{PeriodicOperator}  # (i ad_K)^j/(j+1)! on Kdot, flat over (n,j)
    Heff        :: Vector{QAdd}              # Heff[n+1] = Heff^(n),   n = 0..N-1
    order       :: Int
    wd          :: Num                       # for output reattachment
    t           :: Num                       # for K(t) reconstruction
end
```

The two triangle families are graded pieces of the same operation applied to two different SEEDS,
which is the only thing that distinguishes them. "Dressed" is again the spec's word for
`e^{iK} X e^{-iK}` (eq:kickedA and the `A^(k)_m` recursion), so

    R = e^{iK} H e^{-iK}  -  i e^{iK} d_t e^{-iK}
        \_____________/       \________________/
           dressedH                dressedKdot          (flow-averaged; see the 1/(j+1) weight)

and `R^(n) = sum_j dressedH - sum_j dressedKdot` reads straight off it. They were `A` and `B`,
which carried nothing and left the minus sign looking arbitrary.

The payoff beyond legibility: D6's dressed coupling operator is visibly the SAME routine with the
seed changed to `A`, giving `dressedA`. That is the design's "the seed stays a parameter" claim,
made checkable by the name.

Flat `Vector` rather than `Dict{Tuple{Int,Int}}`, indexed `idx(n,j) = n*(n+1)/2 + j + 1`; see §7 L1
for why. `Heff` is offset by one, since order 0 exists and Julia indexes from 1. That offset is the
one indexing trap here, so it is written out in the field comment and asserted in the accessor.

Stage n (n = 0 .. N-1), given K^(1..n), Kdot^(1..n):

    dressedH^(n)_[0]    = (n == 0) ? H : empty
    dressedH^(n)_[j]    = (i/j)     * sum_{k=1}^{n-j+1} [K^(k), dressedH^(n-k)_[j-1]]
    dressedKdot^(n)_[j] = (i/(j+1)) * sum_{k=1}^{n-j+1} [K^(k), dressedKdot^(n-k)_[j-1]]   j >= 1

    R^(n)      = sum_{j=0..n} dressedH^(n)_[j] - sum_{j=1..n} dressedKdot^(n)_[j]
    Heff^(n)   = time_average(R^(n))
    K^(n+1)    = antiderivative(R^(n) - Heff^(n), gauge)   # mean-subtracted, `DERIVATION.md` §1
    Kdot^(n+1) = derivative(K^(n+1))

`dressedKdot^(n)_[0]` is `Kdot^(n+1)`, not computable at stage n. It is never read:
`dressedKdot^(n)_[j>=1]` reaches only `dressedKdot^(n-k)_[j-1]` with k >= 1, so the order index
strictly decreases. Do not store the (n,0) slot; read `Kdot[n-k+1]` directly, so a demand-driven
memo cannot trip over it.

Key convention: `Kdot[n+1]` holds a *reduced order-n* object. Harmless internally (reduced
`derivative` and `antiderivative` are order-independent) but it is the one place key and order
disagree, and it is where wd gets reattached.

Both families vanish for j > n; the table is O(n^2) nodes. They share one peeling routine, differing
only in seed and in 1/j vs 1/(j+1): `dressedKdot` is `dressedH` averaged over the flow parameter u,
since `e^{iK} d_t e^{-iK} = -i int_0^1 du e^{iuK} Kdot e^{-iuK}` and `int_0^1 u^j du = 1/(j+1)`.

The seed stays a parameter. Two consequences: the arXiv supplement B III graded-H generalization
is a one-line change, and so is the dressed coupling operator (`DECISIONS.md` D6).

## 5. Public API

```julia
struct VanVleck end

vv = floquet_expansion(H::PeriodicOperator, VanVleck(), N)

effective_hamiltonian(vv)        # Heff^[N] = sum_{k<N}, wd reattached
effective_hamiltonian(vv, n)     # order-n piece alone
kick_operator(vv)                # K^[N] = sum_{k<N} K^(k)  :: PeriodicOperator
kick_operator(vv, t)             # K(t) rebuilt with expim phases
```

There is no `kick` option. It was proposed and removed: the triangle builds `dressedH^(n)_[j]` from
`ad_{K^(k)} dressedH^(n-k)_[j-1]`, so every `K^(k)` is required internally whether or not the caller
wants it back. `kick=false` would have saved only the final materialization, i.e. nothing, while reading
as a speedup. Removing it also drops a type parameter and a dispatch-level error path.

The asymmetry is what settles it: adding a keyword with a default later is non-breaking, removing
one is not. And a genuinely K-free engine would be a different GAUGE, hence a different method
type, for which `kick_operator` simply would not be defined. Nothing is foreclosed.

#### No `simplify` kwarg

Proposed and removed. It would have meant "run SQA `simplify` once on the output for readability",
which is cosmetic, while the operation that actually matters is the in-loop prune sweep (§7), which
is about correctness of pruning and about runtime. One boolean must not govern both: wiring
`simplify=false` to the sweep would mean the "fast" setting disables the pruning that makes it
fast, and wiring it only to the output makes it a knob for tidiness that most callers never touch.

Always simplify the output. The sweep is an INTERNAL policy (default: once per order stage). Expose
either as a kwarg only if measurement shows it is wanted; adding a defaulted keyword later is
non-breaking.

`kick_operator(vv)` must EXCLUDE K^(N). Stage N-1 computes K^(N) and Kdot^(N), but the spec's approximant
is `X^[N] = sum_{k<N}`. Easy off-by-one in the accessor, and wasted work at the last stage.

Truncation follows the spec: N=1 is the RWA (H_0 alone), error O(wd^-N). The spec flags that the
Liouvillian literature shifts this by one. See D2.

## 6. Cost

Support: K^(k) has support kM, so order n needs |l| <= (n+1)M. For a finite Fourier polynomial
drive the expansion is exact and finite with no truncation in m.

Dominant cost is SQA commutators; `*` on two QAdds distributes over all |a|x|b| term pairs.

Call NEITHER `normal_order` NOR `simplify` in the loop. SQA's own docstring: `normal_order` is
"the identity on anything built through public arithmetic: *, +, -, ^, commutator", and the
recursion builds everything that way, so calling it re-streams every term for nothing.
`simplify` is documented as a finalizer.

But: the pruning invariant needs care. Operator-level zeros vanish for free (`_addto_key!` deletes
zero coefficients eagerly). Coefficient-level cancellation on the `Complex{Num}` tier does NOT
happen without `simplify`, so mathematically-zero buckets survive as structurally nonzero,
defeating both pruning and the (n+1)M support bound in practice. Needs an explicit periodic
simplify-and-prune sweep with a documented cost knob.

The bucket recursion is the dynamic program collapsing the n-fold index sums of the closed-form
approach into n sequential convolutions (`DERIVATION.md` §6).

## 7. Performance, memoization, inference

Measured, not assumed. SQA 0.10.1, KPO-scale buckets (~4 terms each), 15 harmonics:

    one commutator                    21.9 us
    convolution, Dict-backed  (14x15)  4.32 ms
    convolution, Vector-backed(15x15)  4.36 ms
    210 commutators alone              4.60 ms

**Container choice is irrelevant to speed.** Commutators are ~100% of the time; Dict and dense
Vector are within noise of each other and of the commutator total. So pick the container on
inference and predictability, not throughput, and spend optimization effort on reducing the
NUMBER of commutator calls instead.

### Inference: the repo gates are achievable, contrary to the earlier review

Verified with `Base.return_types` and `@inferred`:

    commutator(::QAdd, ::QAdd)            -> QAdd                        concrete
    (::Rational{Int}) * (::QAdd)          -> QAdd                        concrete
    (::Complex{Rational{Int}}) * (::QAdd) -> QAdd                        concrete
    convolution, Dict-backed              -> Dict{Int,QAdd}   @inferred  OK
    convolution, Vector-backed            -> Tuple{Vector{QAdd},Int}     @inferred OK

`QAdd` is a concrete type; the `Union` lives inside `Coeff.tail`, one level below anything we
return. So CLAUDE.md's `@inferred` gate applies to the whole engine, not just scaffolding. Do not
scope it down. `@allocations` is the one gate that must stay scoped to bookkeeping: building an
SQA expression allocates by construction.

Also verified: `(1//3) * (a'*a)` stays exactly `1//3`, so D5 works as intended.

### Memoization, three levels

**L1, triangle nodes.** Bounds are known up front (n < N, j <= n), so use a flat
`Vector{PeriodicOperator}` with `idx(n,j) = n*(n+1)/2 + j + 1`, not `Dict{Tuple{Int,Int}}`.
Concrete, contiguous, one allocation, no hashing, index-computable. Not faster in any measurable
way, but it is statically sized and trivially inferrable, which the Dict is not.

**L2, freeing intermediates. Not needed; do not build it.** `dressedH^(n)_[j]` reads
`dressedH^(n-k)_[j-1]` for every k >= 1, so every lower-order node stays live to the end. Memory is
O(n^3 * M) QAdds: at n=5, M=8 that is 15 nodes x ~48 buckets = ~720 QAdds; at n=10, ~5000. Both
trivial. A free-as-you-go scheme would buy nothing and risk correctness.

**L3, commutator subexpression memo. Speculative; measure first.** `Dict{Tuple{QAdd,QAdd},QAdd}`
is sound (QAdd hash/isequal are structural and order-independent, and QTerm caches its hash), but
the triangle calls `commutator(K^(k), .)` with a fixed first argument and mostly distinct second
arguments, so pair recurrence is not systematic. Build only if profiling shows repeats. Never
persist such a cache: SQA interning is session-local.

### The optimization that actually pays: Hermiticity halving

Every triangle node is Hermitian. `i*ad_K` maps Hermitian to Hermitian, and the Lie weights are
`(1/j!)(i*ad_K)^j`, so each `dressedH^(n)_[j]` and `dressedKdot^(n)_[j]` is a Hermitian periodic
operator, hence

    X_l^dag = X_{-l}

So compute buckets `l >= 0` only and fill `l < 0` by adjoint. That is a ~2x cut in commutator
calls, which is ~100% of runtime, and it costs nothing. It also doubles as a correctness check:
compute a few negative buckets directly in debug mode and assert they match.

Cheaper still, and already planned: skip zero buckets before calling `commutator` at all.

### Staticness

Do NOT reach for `Val{N}` or `@generated` for the triangle. The loop bounds are runtime `N`, and
specializing per order would multiply compile time for zero runtime gain, since the benchmark shows
loop overhead is unmeasurable next to `commutator`. Keep `N` an ordinary `Int`. Sink type
parameters (`where {F}`) on any helper taking a closure, per CLAUDE.md.

### Invalidations

- No type piracy. Methods are defined only on our own types (`PeriodicOperator`,
  `FloquetExpansion`, `VanVleck`); `Base` methods always have one of ours as an argument.
- Do NOT define scalar-times-operator methods. SQA already covers `Rational`, `Complex{Rational}`
  and `BasicSymbolic` times `QAdd` (verified above). Adding our own would be piracy on two foreign
  types and a genuine invalidation source.
- Precompilation: a `@compile_workload` running an order-2 expansion on a one-mode, two-harmonic
  drive. Per CLAUDE.md this trades precompile time for TTFX; keep the workload small, since the
  runtime cost is all in `commutator`, which precompilation cannot avoid.

### Arithmetic: SETTLED, `Rational{Int}` with a loud ceiling

MEASURED on SQA 0.10.1: coefficient arithmetic does NOT promote on overflow, it throws.

    (1//3037000493)^2 * (a'a)  ->  1//9223371994482243049      exact, at the Int64 ceiling
    (1//3037000493)^3 * (a'a)  ->  OverflowError               9223371994482243049 * 3037000493

So the failure mode is loud, which is what matters: no silent wraparound, no silent precision loss.
Use `Rational{Int}` and do not pre-emptively widen. Also measured: a float-representable rational
collapses to the `Native` (ComplexF64) tier, `1//2^40` printing as `9.09e-13`, but that value is
exactly representable so nothing is lost; non-representable ones such as `1//3` stay exact on the
`Complex{Num}` tier.

Ceiling estimate: denominators carry Deprit `1/j`, `1/(j+1)` and products of harmonic indices up to
`(n+1)M`, so worst case ~`((n+1)M)^n * n!`, which at n=8, M=8 is ~3e19 against the 9.2e18 limit.
That is a WORST case on unreduced fractions; real denominators reduce hard, so the practical
ceiling is higher and must be measured, not predicted. Add a high-order test to find it. If it ever
binds, `Rational{Int128}` at the boundary is the fix, but it costs nothing to defer given the
`OverflowError` is loud.

## 8. Tests

### Algebra layer (`PeriodicOperator`), IMPLEMENTED

Anchored on the one thing the type claims: that it stores `X(t) = sum_l X_l e^{-i l wd t}`. Each
test reconstructs that time-dependent operator and checks the operation against Symbolics' own
calculus, instead of restating the harmonic bookkeeping. Every one was verified to FAIL under a
corrupted implementation before being kept:

| Test | Corruption it was shown to catch |
|---|---|
| `evaluate([K,X]) == [evaluate(K), evaluate(X)]` | bucket index `p-q` for `p+q`; overwrite for sum |
| Jacobi identity on three periodic operators | bucket bookkeeping generally |
| `d/dt evaluate(X) == wd * evaluate(derivative(X))` | sign flip in `derivative` |
| `d/dt evaluate(antiderivative(X)) == wd * evaluate(X)` | sign flip in `antiderivative` |
| `<antiderivative(X, VanVleck())> == 0` | the gauge condition itself |
| `antiderivative` throws on a DC input | the missing mean subtraction (a bug that shipped once) |
| Hermiticity closed under all four ops | precondition of the §7 halving optimization |
| support adds under commutator | the `(n+1)M` bound |
| `49*K[49] == i*a` after `antiderivative` | D5 exactness (float leaves -1.11e-16 at l=49) |
| `@inferred` on all four plus arithmetic | the CLAUDE.md gate |

`d/dt` here is Symbolics `Differential(t)` applied to each coefficient through `to_num`, which
differentiates `expim` correctly. That is what makes the sign tests independent of the convention
they are testing, rather than a restatement of it.

### Collector (`harmonics`), IMPLEMENTED

Anchored on the round trip `harmonics(H, w, t)(w, t) == exponential_form(H)`, in both directions.

| Test | What it holds down |
|---|---|
| round trip on a drive with all of the below at once | the parser as a whole |
| `harmonics(X(w,t), w, t) == X` | the reverse direction, so neither is quietly lossy |
| two harmonics sharing one monomial land in two buckets | F1, the finding that broke the naive design |
| `cos` and `sin` land on +-1 with even/odd symmetry | Base expands a written `exp(im*w*t)` before SQA sees it |
| `cos(w*t + phi)` keeps phi on the coefficient | phi must not be read as part of the harmonic index |
| `w*t^2` and a half-harmonic both throw | accepting either yields a wrong expansion, not an error |
| `3*X[0] == a'a` and `7*X[1] == a` | D5, the parser must not be where exactness is lost |
| a Hermitian drive gives `ishermitian` | eq:fourierH |

PRIMARY GATE for the ENGINE — residual scaling. Needs no oracle, no license, no transcription, and validates K
and Heff jointly at every order:

    finite-dimensional numeric backend, random d x d harmonics with H_{-m} = H_m^dag
    measure  max_t || e^{iK}H_S e^{-iK} - i e^{iK} d_t e^{-iK} - Heff^[N] ||
    across wd in {20,40,80,160}; assert fitted log-log slope = -(N-1) +/- 0.1 for N up to 5-6

The -(N-1) rather than -N is correct: eq:defining0 contains d_t, which costs one order, and the
omitted Kdot^(N) is an order-(N-1) object. Demonstrated sharp against three corruption classes:
(i/(j+1)) -> (i/j), k-range n-j+1 -> n-j, and a sign flip in `antiderivative` flatten every slope
to -1.
Frechet derivative via the [[A,E],[0,A]] block-exponential trick, no extra dependency.

| Level | Target | Source |
|---|---|---|
| unit | eq:K1, eq:R1, eq:Heff1, eq:K2, eq:Heff2 | spec, hand-verified in R3 |
| unit | Heff^(3) (8 terms), Heff^(4) (31 terms) | arXiv:2108.02861 §D transcribed, times (-1)^n |
| e2e | KPO: eq:kpo-heff0, eq:kpo-heff1 | spec example |
| e2e | quantum Kapitza through order 4 | arXiv:2108.02861 B.10 / B.11 |
| unit | QuasienergyOperator eigenvalues (mod wd) vs Heff^[N] spectrum, large wd | nothing else covers it |
| unit | collector round-trip on a hand-built drive | highest-risk component |
| unit | support bound |l| <= (n+1)M | |

DROPPED as vacuous: Hermiticity of Heff (S and W are each separately Hermitian, so any real
coefficients pass; confirmed to hold even for corrupted recursions) and `<K^(n)> == 0` (true by
construction, since `antiderivative` under the van Vleck gauge drops bucket 0). Neither can catch
a wrong weight.

Comparing against the published tables needs a harness: their entries are index-sums over free m_i
with constraints (m_i != 0, m_1 != m_2, ...), ours are concrete integer buckets. Instantiate their
formula on a specific H and expand the sums, transcribing constraints exactly. That harness is
where errors hide, which is another reason the residual test is the primary gate.

Oracle translation. Define T = [substitute H^V_m -> H_{-m}, then flip every summation index
m_i -> -m_i, then hbar -> 1]. Then, uniformly:

    Heff^(n) =  T( their K^(n) )
    K^(n)    = -T( their S^(n) )

Both verified at n=1 and n=2. RESOLVED: the kick rule is a constant -1 once T is applied. Read
against the printed form directly it looks graded, `(-1)^(n+1)`, because T itself contributes
`(-1)^n` (denominators are homogeneous of degree n and flip sign; numerators return to printed form
since all index arguments are linear homogeneous). State it as T to avoid the confusion.

Known paper defect: their S^(3) fifth term prints `3 m1^3`; correct is `3 m1^2 m2` (as printed it
vanishes identically).

Licensing: github.com/xiaoxuisaac/vanVleck-recursion has NO LICENSE, so all rights reserved.
Oracle values in tests come from the published paper (citable transcription) only. Regenerating
them from our own engine would be a regression test, not a correctness test, so it does not
substitute. This is why the residual test carries the load.

Comparison rule is `iszero(simplify(a - b))`. Note this requires exact weights (`DECISIONS.md` D5): it cannot
match the spec's exact rationals in eq:kpo-heff1 (17/4, 3/20, 9/40) if weights are floats.

Repo gates (`@inferred`, `@allocations`, JET, ExplicitImports): scope them to the bookkeeping
scaffolding (key construction, bucket iteration, the (n,j) loop), not to anything returning an SQA
expression, which allocates by construction. A blanket application will just get disabled.

Re-export note: `expim`, `exponential_form`, `trigonometric_form` are `@public`, NOT exported, so
`@reexport using SecondQuantizedAlgebra` does not forward them. Users need them for input
normalization and for reading `kick_operator(vv, t)`. Re-export explicitly.

## 9. Cheap structural checksums

Before comparing values against the published tables, compare structure. From R2, term counts:

    order      0    1    2    3    4    5
    Heff       1    1    2    8   31  138
    generator  -    1    2    6   21   84

A term-count mismatch localises an error to an order in one cheap assertion, before any coefficient
comparison. R2 also recorded coefficient-sum checksums at order 5 (Heff: -4583/720 in their sign
convention; generator: -11) and at lower orders (Heff^(3): -37/24, Heff^(4): 367/120), which catch
sign and weight errors without term-by-term matching.

Transcription warning carried over from R2: in the printed order-5 generator block, negative terms
appear as a double `+ -` glyph, not a single `-`. It means negative. Do not read it as `+/-`.

## 10. The series is asymptotic, and users must be told

Giorgilli and Sansottera measure `||Phi_s|| ~ s!` for normal-form remainders, confirmed numerically
by root and ratio criteria. This is intrinsic to the method, not an implementation defect. The
practical consequence is that **optimal truncation order is finite and problem-dependent**: for one
and the same Hamiltonian they report optimal order 32 at moderate energy and order 9 at high
energy.

A package advertising "arbitrary order" without saying this is misleading. Two obligations:
- Docs must state the expansion is asymptotic, and that going to higher order past the optimal
  truncation makes the answer WORSE, not better.
- Consider exposing a per-order norm or term-count so a user can see the series turning. The
  residual-scaling test machinery (§8) already computes exactly the quantity needed.

The spec's own truncation paragraph says what the truncation costs "is priced separately", so this
is consistent with, not additional to, the source.

## 11. Risks

**R-1, collector index extraction. RETIRED, with no residual.** Prototyped and working; see §3.
Extraction runs through `to_num` plus a SymbolicUtils walk, with no `Poly`/`Monomial` access, and
term iteration is the public Base protocol on `QAdd` (F4). Parsing happens once on the seed drive,
never in the loop. Nothing here depends on an SQA internal, and no upstream change is required.

**R-2, term growth.** Nested commutators compound SQA's O(|a|x|b|) product. Benchmark at order 3
against the KPO early. Fallback is memoizing commutator subexpressions on QAdd (safe: hash/isequal
are structural and order-independent; never persist, since interning is session-local), capping m,
or caching per-H results. NOT precomputed universal-formula tables: that is the approach `DECISIONS.md` argues
against, it carries (m1-m2) denominators the bucket design never produces, and sourcing it from the
reference implementation would violate the licensing rule.


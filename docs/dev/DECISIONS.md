# Decisions and rejected paths

Why the design is what it is. The design itself is in `DESIGN.md`.

Read this before re-proposing anything: several reasonable-looking options were evaluated and
rejected, and the reasons are recorded so they are not re-litigated.

## Decision log


D1. Collector takes BOTH a parsed `QAdd` and a direct Dict. Parsed path documented as primary;
    the direct constructor keeps the engine testable without the parser.

D2. Spec convention: `X^[N] = sum_{k<N}`. `floquet_expansion(H, VanVleck(), 1)` is the RWA.
    Document loudly; the Liouvillian literature counts shifted by one.

D3. wd reattached into a single expression by default, per-order access via
    `effective_hamiltonian(vv, n)`.

D4. Names and method options:
      VanVleck()                           plain singleton; no `kick` option (proposed, removed)
      PeriodicOperator                     the harmonic container (Toeplitz/periodic subclass)
      QuasienergyOperator(H, nmax)         Q = H_S - i d_t, a real type, NOT Toeplitz
      floquet_expansion(H, VanVleck(), N)  entry point, SciML `solve(prob, alg, args...)` order
      FloquetExpansion{VanVleck}           result; structure shared with Floquet-Magnus,
                                           only the gauge fixing differs
      effective_hamiltonian / kick_operator / harmonics
    No `van_vleck` alias. `expand` was rejected: it collides with SQA's exported `expand`.

    Internal names come from the SPEC's prose, not invented. Checked against sec:vv and sec:comb:
      antiderivative(X, gauge)  was `antiderivative0`. The spec says "zero-average antiderivative"
                                (sec:vv, on eq:recursion-split), so `antiderivative` is its word;
                                the `0` was the "zero-average" qualifier, i.e. the van Vleck gauge
                                encoded as a digit, which becomes a lie under Floquet-Magnus (D7).
                                Promote the qualifier to an argument. Rejected: `solve_homological`
                                (the spec does call this the homological equation, but the name does
                                not say which direction it is solved), `integrate` (collides with
                                `average`, itself a time integral), `primitive` (collides with the
                                design's own "four primitives").
      time_average(X)           was `average`. SQA EXPORTS `average` for the expectation value
                                (`average(a'*a)` is `⟨a' * a⟩`), and we `@reexport` SQA, so that
                                name is in every user's scope with a different meaning. Contrast
                                `commutator`, which IS the same operation on a new type and so
                                extends `SQA.commutator`.
      dressedH, dressedKdot     were the triangle families `A` and `B`. "Dressed" is the spec's word
                                for `e^{iK} X e^{-iK}` (eq:kickedA), so the two families are named
                                by their seed and `R = sum dressedH - sum dressedKdot` maps onto
                                `R = e^{iK} H e^{-iK} - i e^{iK} d_t e^{-iK}` term by term. With
                                A/B the minus sign in R looked arbitrary. Bonus: D6's dressed
                                coupling operator is then `dressedA`, the same routine reseeded,
                                which is exactly what D6 claims. Rejected: `Hrot`/`Kdotrot` ("rot"
                                is curl to half the audience, and neither is a rotation).

D5b. Coefficient type SETTLED as `Rational{Int}` (measured, SQA 0.10.1). Overflow THROWS
    `OverflowError` rather than promoting or wrapping, so the ceiling (~order 8 worst case) is a
    loud failure, not a silent corruption. Do not pre-emptively widen to Int128/BigInt; measure
    where it actually binds first, since real denominators reduce far below the worst case.

D5. EXACT rationals for the Deprit weights. Accepts the slow coefficient tier from order 2 up.
    Required for `iszero` reliability and to reproduce eq:kpo-heff1's 17/4, 3/20, 9/40.

    MEASURED NUANCE: exactness is preserved but DISPLAY is not uniform. SQA's `Native` tier is
    ComplexF64, so a float-representable rational collapses to it: `3//4` prints as `0.75`, while
    `1//3` stays on the `Complex{Num}` tier and prints as `1//3`. Both are exact, and
    `iszero(simplify(a-b))` is unaffected, but eq:kpo-heff1's `17/4` will print as `4.25` while
    `3/20` prints as `3//20`. Decide whether to post-process the output into uniform rational
    display; it is presentation only, not correctness.

D7. v1 is scoped to van Vleck. The method type already names the choice, so later additions need
    no restructuring. But they sit on TWO DIFFERENT AXES and must not be conflated:

    AXIS 1, gauge within the same factorization  U = e^{-iK(t)} e^{-i Heff (t-t')} e^{+iK(t')}
      VanVleck        <K> = 0          bucket 0 dropped                v1
      FloquetMagnus   K(t_0) = 0       constant fixed by evaluating at t_0
    Differ ONLY in the integration constant of `antiderivative`. Same triangle, same collector,
    same everything else. One line. Make the constant a property of the gauge type NOW so this is
    an added method rather than a refactor.
    Note FloquetMagnus's K has NONZERO mean, so its generator is not block-off-diagonal in Sambe
    space. It is a gauge of the factorization, not of a block diagonalization.

    AXIS 2, a different perturbation SCHEME
      DesCloizeaux    Bloch wave operator + symmetric orthonormalization,
                      H_dC = (P Om^dag Om P)^{-1/2} P Om^dag H Om P (P Om^dag Om P)^{-1/2}
    NOT reachable by changing an integration constant. Different recursion, different algorithm.
    Shares only the `PeriodicOperator` layer and the collector. Budget it as an engine, not a flag.
    It enters the plan only because Pymablock may turn out to be this rather than van Vleck
    (see "Rejected" below); if so it ships here rather than being discarded.

D6. Dressed coupling operator OUT OF SCOPE. The seed stays parameterized regardless, since that
    costs nothing and is what makes both this and the graded-H case one-line additions later.

## Rejected: Pymablock and Nikolaev as a faster engine


Pymablock's recursion (arXiv:2404.03728, SciPost Phys. Codebases 50). **BSD-2-Clause**, so unlike
the unlicensed vanVleck Python code we may read the source and port the algorithm with attribution.
That difference matters: `DESIGN.md` §8's licensing rule does NOT apply here.

#### The Sambe translation (derived here; this is the part their paper does not hand us)

Their split is unperturbed `H_0` plus perturbation `H'`. Ours, in Sambe space:

    Q[m,n] = H_{m-n} - m*wd*delta_mn
             \_______/  \___________/
             perturbation   H_0

The `-m*wd` diagonal is the LARGE part, since the expansion parameter is 1/wd. The Toeplitz part
`H_{m-n}` is the perturbation. "Block-diagonal" means `m = n`, i.e. harmonic 0 only, which is
exactly `<.>`.

Sylvester: `E_m = -m*wd`, so `E_m - E_n = -(m-n)*wd = -l*wd`. The solve is a division by `-l*wd`,
i.e. our `antiderivative` up to the factor of i that comes from `-i d_t` versus `d_t`.

**Everything stays a PeriodicOperator; we never materialize a Sambe matrix.** Proof: for Toeplitz X,

    ([H_0, X])[m,n] = -m*wd*X_{m-n} + n*wd*X_{m-n} = -(m-n)*wd * X_{m-n}

which is Toeplitz, and products of Toeplitz operators are Toeplitz (that is the convolution we
already have). So the recursion closes inside `PeriodicOperator` and reuses `DESIGN.md` §2's four operations
plus one new primitive, the Sylvester division. **No new container, no new algebra.**

#### The recursion

`U = 1 + U'` with `U' = W + V`, W Hermitian, V anti-Hermitian. Unitarity alone fixes W:

    U^dag U = 1 + U' + U'^dag + U'^dag U' = 1
    U' + U'^dag = 2W                              (V is anti-Hermitian, cancels)
    => W = -U'^dag U' / 2

This is a legal recurrence because `U'` has no 0th-order term, so `W` at order n reads `U'` only at
orders < n. **One Cauchy product per order, hence O(n) products per order and O(n^2) total.**

V is fixed per order by requiring the l != 0 part of the transformed operator to vanish, which is
one Sylvester solve, free here.

Their design constraints, worth stating as invariants to assert:
- exactly one Cauchy product per order;
- **never multiply by H_0** — their words: *"any additional multiplications by H_0 must cancel with
  additional energy denominators… unnecessary work, and it gives longer intermediate expressions"*.
  In our terms: never multiply by the `-l*wd` diagonal. If an implementation does, it is wrong;
- exactly one Cauchy product by the selected (harmonic-0) part.

#### TODO before implementing

The auxiliary series `X = [U', .]` and the CSE'd final form (their Eq. 31, with `A = H'_R U'` and
`B = X - H'_R - A`, and `U'^dag X` computed once and reused via `+h.c.`) must be **transcribed from
the paper and their source**, not reconstructed from a summary. Their CSE is worth the trouble:
orders 2/3/4 cost 1/3/11 products versus 1/4/27 for an optimally-implemented reference. Read
`pymablock/block_diagonalization.py` and `pymablock/algorithms.py`.

Also read `pymablock/number_ordered_form.py` before finalizing anything about our canonical form.
They key terms by signed powers per mode and allow values containing non-polynomial functions of
number operators such as `1/(N+1)`, which a normal-ordered form cannot represent at all. We are
committed to SQA's normal-ordered form, but we should know what that costs us.

Their `BlockSeries` is a lazy memoized object with masked structural zeros, and the optimizations
are emitted by a DSL plus AST code generator because they would otherwise swamp the algorithm's
definition. Signal: **the complexity lives in the bookkeeping, not the math.** If our version starts
to sprawl the Julia analogue is a macro over a declarative spec, not hand-written special cases.

#### Why we do NOT port this as the fast van Vleck path

Question raised and worked through: if the difference is "just a gauge", can we run Pymablock's
O(n^2) recursion in the van Vleck gauge? No, and the reason matters.

The cost is not in the gauge, it is in **what the recursion is parameterized by**:

    Pymablock   condition lands on U itself   (V, the anti-Hermitian part of U-1, is
                                               block-off-diagonal), so W follows free from
                                               unitarity as one Cauchy product
    van Vleck   condition lands on log U      (the GENERATOR is block-off-diagonal)

Those agree at leading order and diverge later, which is exactly the order-3 split measured below.
Converting between them is the Theta(n^3) log. So imposing van Vleck on their recursion dismantles
the thing that makes it fast; it is not a parameter change.

**The decisive point is the generator, not the asymptotics.** The current standings, after reading
Nikolaev in full (an earlier draft claimed Nikolaev "fixes the triangle" to O(n^2),
which was a misreading and is corrected there):

    algorithm              total cost   gauge              gives the generator?
    Deprit triangle (ours) O(N^3)       van Vleck, VERIFIED    yes
    Nikolaev 6.1 + 6.2     O(N^2)       plausibly van Vleck    yes, directly
    Pymablock              O(N^2)       unverified, at risk    NO

**DECISION: Pymablock is not the fast van Vleck path, because it does not produce the generator at
all** and recovering it via log costs Theta(n^3), which is the triangle's cost anyway. If a faster
van Vleck exists it is Nikolaev's explicit generator, which produces K directly and whose gauge
argument is much stronger (its reduced resolvent has zero kernel component by construction, i.e.
`<W> = 0`).

Speed work that is safe regardless of which engine wins:
  - Hermiticity halving                      ~2x on ~100% of runtime
  - zero-bucket pruning before the commutator call
  - Pymablock's CSE *ideas* where they transfer (subexpression reuse is gauge independent)

Pymablock's remaining value is narrow and none of it is on the v1 path: CSE technique to borrow,
a cross-check if its gauge turns out to match, and a candidate gauge to ship later.

#### Scope: v1 is van Vleck gauge only

Decided. v1 ships the van Vleck gauge; other gauges are a later iteration. Two consequences:

1. **The gauge question no longer blocks v1.** The triangle gives the van Vleck gauge by
   construction, since `int0` drops bucket 0 so `<K> = 0` identically. It only gates the fast path.
2. **If Pymablock turns out to be des Cloizeaux, nothing is wasted.** It becomes `DesCloizeaux()`
   in the later gauge iteration rather than a discarded engine. Build it either way; only the name
   it ships under is contingent.

So: v1 = triangle only, under `VanVleck`. A Pymablock engine ships later as its own gauge type if
the check below says its gauge is not van Vleck's, or as a cross-check if it is.

#### Gauge prerequisite — TESTED, and the risk is REAL

Attempted to settle this numerically and got a negative result worth having.

Test: van Vleck via the Deprit triangle versus an exact canonical block diagonalization of the
truncated Sambe operator, using the des Cloizeaux direct rotation
`U = (Pt P + (1-Pt)(1-P)) [X^dag X]^{-1/2}`. Random Hermitian d=3 harmonics, |m| <= 2, comparing
operators (not spectra: spectra are gauge invariant and cannot discriminate).

    N=1  slope -1.000   expect -1     agree
    N=2  slope -2.005   expect -2     agree
    N=3  slope -2.986   expect -3     agree
    N=4  slope -3.020   expect -4     DISAGREE
    N=5  slope -3.017   expect -5     DISAGREE

The floor is NOT truncation: identical at Mcut = 14, 24, 40, and N=4 and N=5 give the same residual.
So van Vleck and des Cloizeaux agree through order 2 and **diverge at order 3**.

Diagnostics on the des Cloizeaux construction:
- Heff is translation covariant, `||Heff(blk 0) - Heff(blk 1)|| ~ 1e-13`, as Floquet symmetry
  requires.
- But its generator is **not Toeplitz**: block deviations 4.0e-2 at |l|=1 and 3.1e-2 at |l|=2
  against `||G|| = 5.8e-2`. It is the transformation tailored to one target block, not the global
  periodic K(t).

This matches established theory: van Vleck, des Cloizeaux and Bloch are DISTINCT quasi-degenerate
perturbation schemes. van Vleck and des Cloizeaux both give Hermitian effective Hamiltonians and
agree at low order, but they are not the same transformation.

**Consequence: "canonical block diagonalization" is not a single thing, and the schemes diverge at
exactly order 3.** The prerequisite is therefore a real correctness risk, not a formality, and it
cannot be discharged ahead of transcribing Pymablock's specific construction. Their Floquet
formulation uses a bosonic ladder (`Omega*N_a`, couplings via `a`, `a^dag`), which is built from
translation-covariant objects and so is PLAUSIBLY the Toeplitz/van Vleck gauge rather than the
des Cloizeaux one, but plausibly is exactly the standard this test just failed.

Discharge it like this, in order:
1. transcribe Eq. 31 and their `block_diagonalization.py`;
2. run their construction against eq:Heff2 (order 2 will NOT discriminate, per the table above)
   and against Heff^(3) from `DESIGN.md` §9's checksums, which will;
3. it ships under its own method type either way. If its gauge is van Vleck's it is a cross-check
   and a source of CSE technique; if it is des Cloizeaux-like it is a later gauge offering. Its
   Heff must never be presented as van Vleck's without the order-3 check passing.

Both diagnostics above should become standing tests: the order-scaling comparison, and the
generator-Toeplitz check (`G[m,n]` must depend only on `m-n` for a genuine van Vleck generator).

Note that des Cloizeaux is itself fully symbolic (Bloch wave operator plus symmetric
orthonormalization; see RESEARCH.md §5), so it is a rival ENGINE, not just a numerical yardstick.
Pymablock's `W = -U'^dag U'/2` carries the symmetric-orthonormalization signature. Treat the gauge
risk as likely rather than remote.

An attempt to shortcut the question by Toeplitz-projecting the des Cloizeaux generator FAILED and
is recorded in RESEARCH.md §5 so nobody repeats it. Separate two conditions when retrying:
(a) block-off-diagonal, the SW/van Vleck condition, and (b) Toeplitz, the Floquet condition.
van Vleck needs both; the experiment measured only (b).

#### What the two engines share

The collector, `PeriodicOperator` and its four operations, the wd-reduced convention, exact
rationals, the numeric backend, and every test in `DESIGN.md` §8. Only the recursion differs. Keep the shared
layer free of any engine-specific assumption.

#### Cross-validation

Once both exist and the gauge question is settled, running them against each other is **the
strongest correctness test in the plan**: two structurally unrelated algorithms, sharing no code
path beyond the SQA algebra, computing the same object, at every order, for free. Stronger than any
transcribed oracle. Make it a standing test, not a one-off check.


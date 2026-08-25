# Derivation: arbitrary-order van Vleck recursion, Fourier domain

The mathematics behind `DESIGN.md`. Conventions are the spec's.

Spec: tmp/floquet_bath_tracing.tex sec:vv. All notation below is the spec's, not arXiv:2108.02861's.

## 0. Conventions

    H_S(t) = sum_{m in Z} H_m e^{-i m wd t},   H_{-m} = H_m^dag
    U_S(t,t') = e^{-iK(t)} e^{-i Heff (t-t')} e^{+iK(t')}
    K(t+T) = K(t),  <K> = 0            (van Vleck gauge)
    ad_K X = [K, X]

Order counting: X^(n) = O(wd^-n). K = sum_{k>=1} K^(k). H_S = O(wd^0).
d/dt lowers the order by one: Kdot^(k) = O(wd^-(k-1)).

Defining identity (eq:defining):

    Heff  =!  sum_{j>=0} (i^j / j!)     ad_K^j H_S(t)
            - sum_{j>=0} (i^j / (j+1)!) ad_K^j Kdot(t)

## 1. Residual at arbitrary order

Expand K = sum_k K^(k) and collect wd^-n. The unknown K^(n+1) appears only in the
second sum at j=0 (where k_0 = n+1 is forced), contributing -Kdot^(n+1).
Everything else is the residual:

    R^(n)(t) =   sum_{j>=0} (i^j/j!)     sum_{k_1+...+k_j = n}       ad_{K^(k_1)}...ad_{K^(k_j)} H_S(t)
               - sum_{j>=1} (i^j/(j+1)!) sum_{k_0+k_1+...+k_j = n+1} ad_{K^(k_1)}...ad_{K^(k_j)} Kdot^(k_0)

all indices >= 1. In the second sum j>=1 forces k_0 <= n, so only known kicks enter.

Then (eq:recursion-split):

    Heff^(n) = <R^(n)>
    K^(n+1)  = antiderivative[ R^(n) - <R^(n)>, gauge ]

where `antiderivative` inverts d/dt and the gauge supplies the integration constant. Van Vleck
takes it to be zero, i.e. <K^(n+1)> = 0.

## 2. Deprit triangle (the efficient form)

The composition sums above are exponential in n. Factor them by bracket depth j.
Two families, both memoized on (n, j):

    dressedH^(n)_[j]    := (i^j/j!)     sum_{k_1+..+k_j=n}   ad_{K^(k_1)}..ad_{K^(k_j)} H_S
    dressedKdot^(n)_[j] := (i^j/(j+1)!) sum_{k_0+..+k_j=n+1} ad_{K^(k_1)}..ad_{K^(k_j)} Kdot^(k_0)

Peeling the outermost bracket gives the recursions

    dressedH^(n)_[0]    = delta_{n,0} H_S
    dressedH^(n)_[j]    = (i/j)     sum_{k=1}^{n-j+1} ad_{K^(k)} dressedH^(n-k)_[j-1]

    dressedKdot^(n)_[0] = Kdot^(n+1)             (the unknown; excluded from R)
    dressedKdot^(n)_[j] = (i/(j+1)) sum_{k=1}^{n-j+1} ad_{K^(k)} dressedKdot^(n-k)_[j-1]

    R^(n) = sum_{j=0}^{n} dressedH^(n)_[j] - sum_{j=1}^{n} dressedKdot^(n)_[j]

Both vanish for j > n, so the table is genuinely triangular: O(n^2) nodes.
dressedKdot^(n)_[1] needs dressedKdot^(n-k)_[0] = Kdot^(n-k+1) with k>=1, hence order <= n.
The triangle closes.

This is the (n,k) grid of arXiv:2108.02861 Fig. 2, restated in the spec's gauge and sign conventions.

## 3. Fourier-domain realization

Every object is a harmonic dict {l -> X_l}, X_l a QAdd (SQA 0.10 has no QNumber and no QMul). The four primitives:

    (ad_{K} X)_l = sum_p [K_p, X_{l-p}]          convolution
    (dX/dt)_l    = -i l wd X_l
    <X>          = X_0
    antiderivative(X)_l    = (i/(l wd)) X_l for l != 0,  0 for l = 0   (van Vleck gauge)

Check: d/dt[(i/(l wd)) X_l e^{-i l wd t}] = X_l e^{-i l wd t}. Correct.

Harmonic support: K^(k) has support kM, so ad_{K^(k1)}...ad_{K^(kj)} H_S reaches
(sum k_i)M + M = (n+1)M independent of j. Order n needs |l| <= (n+1)M.
(NOT (j+1)M: that intermediate claim was wrong. The (n+1)M conclusion stands.) Finite and exact for a finite Fourier polynomial drive.
No truncation in m is required in that case; it is only needed for a drive with
unbounded harmonic content.

## 4. Verification against the spec's closed forms

n=0: only j=0 survives both sums. R^(0) = H_S. Then
     Heff^(0) = H_0                                                    [eq:K1] OK
     K^(1)_m  = i H_m/(m wd),  m != 0                                  [eq:K1] OK

n=1: A: j=1,k_1=1 -> i[K^(1),H_S].  B: j=1,k_0=k_1=1 -> -(i/2)[K^(1),Kdot^(1)].
     R^(1) = i[K^(1),H_S] - (i/2)[K^(1),Kdot^(1)]                      [eq:R1] OK

     In Fourier, with Kdot^(1)_m = H_m (m!=0):
     R^(1)_l = -[H_l,H_0]/(l wd) * 1[l!=0]
               - (1/2) sum_{m!=0, m!=l} [H_m, H_{l-m}]/(m wd)         [eq:R1] OK

     Heff^(1) = R^(1)_0 = -(1/2) sum_{m!=0} [H_m,H_{-m}]/(m wd)
                        =  sum_{m>0} [H_m^dag, H_m]/(m wd)             [eq:Heff1] OK
     K^(2)_l  = (i/(l wd)) R^(1)_l
              = -i{ [H_l,H_0]/(l^2 wd^2)
                    + (1/2) sum_{m!=0,l} [H_m,H_{l-m}]/(l m wd^2) }    [eq:K2] OK

n=2: A: j=1,k=2 -> i[K^(2),H_S];  j=2 -> -(1/2)[K^(1),[K^(1),H_S]]
     B: j=1,(k_0,k_1)=(1,2),(2,1) -> -(i/2)([K^(2),Kdot^(1)] + [K^(1),Kdot^(2)])
        j=2 -> +(1/6)[K^(1),[K^(1),Kdot^(1)]]

     Period-average each with <[X,Y]> = sum_l [X_l, Y_{-l}]. Define

       S := sum_{l!=0} [[H_{-l},H_0],H_l]/(l^2 wd^2)
       W := sum_{m!=0} sum_{m'!=0,m} [[H_{-m},H_{m-m'}],H_{m'}]/(m m' wd^2)

     The five contributions collapse (using [X,[Y,Z]] = -[[Y,Z],X] and relabelling)
     to coefficients 1/2 on S and 1/4 + 1/12 = 1/3 on W:

       Heff^(2) = (1/2) S + (1/3) W                                    [eq:Heff2] OK

     The 1/3 is the sharp test: it comes from -1/4 + 1/2 - 1/6 = 1/12 on one
     rearrangement plus 1/4 on the other. Getting the Deprit weights or the
     [X,[Y,Z]] antisymmetry wrong breaks it.

## 5. Consequences for the design

- Memo table keys are structural, (n,j), not operator expressions. Whether
  SecondQuantizedAlgebra can hash expressions is therefore NOT on the critical path.
- Needed algebra on the harmonic series: convolution-commutator, d/dt, <.>, antiderivative.
  Nothing else.
- K^(k) and Kdot^(k) are both wanted at every order; store Kdot alongside K rather
  than recomputing (it is just the diagonal rescaling -i l wd).
- dressedH and dressedKdot share the identical peeling recursion apart from the 1/j vs 1/(j+1)
  prefactor and the seed. One generic routine, two seeds.

## 6. Cross-check against arXiv:2108.02861 Section D (oracle)

Their conventions: H(t) = sum_m H^V_m e^{+i m w t}; K^(n) = static Kamiltonian;
S^(n) = generator; L_S = [S,.]/(i hbar); gauge S_avg = 0 (= van Vleck). Denominators (hbar w)^n.

Translation to the spec's conventions:
  1. same physical H(t) requires  H^V_m = H_{-m}
  2. then flip every summation index m_i -> -m_i
  3. hbar -> 1

Every index argument is linear in the m_i, so step 1 followed by step 2 returns each
H subscript to its printed form. Only the denominator changes sign, and the denominator
is homogeneous of degree n at order n. Hence the whole table translates by one factor:

    Heff^(n) = (-1)^n * [ their K^(n), hbar=1, symbols read as our H_m ]
    K^(n)    = -T( their S^(n) )      RESOLVED, verified at n=1 and n=2

where T is the translation itself: substitute H^V_m -> H_{-m}, flip every summation index
m_i -> -m_i, set hbar = 1. Stated that way both rules are clean:

    Heff^(n) =  T( their K^(n) )
    K^(n)    = -T( their S^(n) )

n=1 check. Printed: S^(1) = i sum_{m1!=0} H^V_{m1} e^{+i m1 w t} / (m1 w).
  substitute -> i sum H_{-m1} e^{+i m1 w t}/(m1 w); relabel m1 = -l
              -> -i sum_l H_l e^{-i l w t}/(l w)  =  -K^(1).            OK

n=2 check, first term. Printed S^(2)/(i) has [H^V_m1, H^V_0]/(m1^2 w^2) e^{+i m1 w t}.
  T gives +i [H_l, H_0]/(l^2 w^2) e^{-i l w t}; eq:K2's first term is -i [H_l,H_0]/(l^2 w^2). OK
  Second term likewise: T gives +i [H_m, H_{l-m}]/(2 l m w^2), eq:K2 has -i times that.  OK

Read against the PRINTED form directly the kick rule looks graded, (-1)^(n+1), because T itself
contributes (-1)^n: denominators are homogeneous of degree n and flip sign, while numerators return
to printed form because every index argument is linear homogeneous. Always state the rule as T to
avoid that confusion.

Verified:
  n=1: (-1)*[H_m1,H_-m1]/(2 m1)  =  -(1/2) sum_m [H_m,H_-m]/(m wd)     matches eq:Heff1
  n=2: (+1)*{ [[H_m1,H_0],H_-m1]/(2 m1^2) + [[H_m2,H_m1-m2],H_-m1]/(3 m1 m2) }
       matches eq:Heff2 under m1 -> -m, m2 -> -m'                       matches eq:Heff2

So their order-3 (8 terms) and order-4 (31 terms) tables are usable as oracles for
Heff^(3), Heff^(4) after multiplying by (-1)^3 and (-1)^4.

Known defect in the paper: S^(3) fifth term prints denominator 3 m1^3; it must be
3 m1^2 m2. As printed the term vanishes identically (sum over m2 of [H_m2,H_-m2] is
odd under m2 -> -m2). Do not use that term as an oracle uncorrected.

Their reference implementation is github.com/xiaoxuisaac/vanVleck-recursion.
NO LICENSE FILE, so all rights reserved: it may be run to generate oracle values,
but no code and no vendored output may be committed to this package.

Their data structure is NOT harmonic-indexed: a term carries only rotating in {0,1}
plus an integration count, and the Fourier indices m_i are synthesised at print time.
That is why their output contains denominators like (m1-m2) and (m1-m2)^2. Our
harmonic-bucket design never produces those: the index is a concrete Int per bucket,
so the denominator is always a plain integer. Different goal (their closed form is
general in H; ours evaluates a specific H), and ours is the one the package wants.

Extension hook (their supplement B III): grading H itself by order, i.e. seeding
dressedH^(n)_[0] = H_S^(n) for all n instead of only n=0. Out of scope now; the seed
should stay a parameter so this is a one-line change later.

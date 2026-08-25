using LinearAlgebra, Random, Printf
const HS = Dict{Int,Matrix{ComplexF64}}

zed(d) = zeros(ComplexF64, d, d)
function addto!(o::HS, l, m); o[l] = haskey(o,l) ? o[l]+m : copy(m); end
function prune(x::HS); HS(l=>v for (l,v) in x if norm(v) > 0); end
function comm(K::HS, X::HS)
    o = HS()
    for (p,kp) in K, (q,xq) in X; addto!(o, p+q, kp*xq - xq*kp); end
    prune(o)
end
scal(c, X::HS) = prune(HS(l=>c*v for (l,v) in X))
plus(A::HS,B::HS) = (o=HS(); for (l,v) in A addto!(o,l,v) end; for (l,v) in B addto!(o,l,v) end; prune(o))
minus(A::HS,B::HS) = plus(A, scal(-1,B))
deriv(X::HS, w) = prune(HS(l=>(-im*l*w)*v for (l,v) in X))
avg(X::HS, d) = get(X, 0, zed(d))
int0(X::HS, w) = prune(HS(l=>(im/(l*w))*v for (l,v) in X if l != 0))

# --- van Vleck via the Deprit triangle (my recursion) ---
function vanvleck(H::HS, N::Int, w, d)
    K = Dict{Int,HS}(); Kd = Dict{Int,HS}(); Heff = Dict{Int,Matrix{ComplexF64}}()
    A = Dict{Tuple{Int,Int},HS}(); B = Dict{Tuple{Int,Int},HS}()
    Bget(m,j) = j == 0 ? get(Kd, m+1, HS()) : get(B, (m,j), HS())
    for n in 0:N-1
        A[(n,0)] = n == 0 ? H : HS()
        for j in 1:n
            acc = HS()
            for k in 1:(n-j+1); acc = plus(acc, comm(get(K,k,HS()), get(A,(n-k,j-1),HS()))); end
            A[(n,j)] = scal(im/j, acc)
        end
        for j in 1:n
            acc = HS()
            for k in 1:(n-j+1); acc = plus(acc, comm(get(K,k,HS()), Bget(n-k, j-1))); end
            B[(n,j)] = scal(im/(j+1), acc)
        end
        R = HS()
        for j in 0:n; R = plus(R, get(A,(n,j),HS())); end
        for j in 1:n; R = minus(R, get(B,(n,j),HS())); end
        Heff[n] = avg(R, d)
        K[n+1]  = int0(minus(R, HS(0=>avg(R,d))), w)
        Kd[n+1] = deriv(K[n+1], w)
    end
    Heff
end

# --- exact block diagonalization of the Sambe operator, canonical (des Cloizeaux) gauge ---
function sambe_heff(H::HS, w, d, Mcut)
    ms = -Mcut:Mcut; nb = length(ms); D = d*nb
    Q = zeros(ComplexF64, D, D)
    blk(i) = ((i-1)*d+1):(i*d)
    for (i,m) in enumerate(ms), (j,n) in enumerate(ms)
        haskey(H, m-n) && (Q[blk(i), blk(j)] .+= H[m-n])
        m == n && (Q[blk(i), blk(j)] .-= m*w*I(d))
    end
    Q = (Q + Q')/2
    i0 = findfirst(==(0), collect(ms)); rows0 = blk(i0)
    P = zeros(ComplexF64, D, D); P[rows0, rows0] = I(d)
    F = eigen(Hermitian(Q))
    ov = [norm(P*F.vectors[:,i]) for i in 1:D]
    sel = sortperm(ov, rev=true)[1:d]
    V = F.vectors[:, sel]; Pt = V*V'
    X = Pt*P + (I-Pt)*(I-P)
    U = X * inv(sqrt(Hermitian(X'X)))          # direct (canonical) rotation
    Heff = (U'*Q*U)[rows0, rows0]
    (Heff + Heff')/2
end

Random.seed!(20260825)
d = 3; M = 2
H = HS()
for m in 0:M
    A = randn(ComplexF64, d, d)
    H[m] = m == 0 ? (A+A')/2 : A
end
for m in 1:M; H[-m] = H[m]'; end


# ---- Nikolaev operators, reduced units (wd = 1) ----
# L_{H0} X |_l = i*l*X_l   =>   S = pseudo-inverse, P = kernel projector (l = 0)
Pk(X::HS, d) = HS(0 => avg(X, d))
Sk(X::HS)    = prune(HS(l => v/(im*l) for (l,v) in X if l != 0))
Lk(F::HS, G::HS) = scal(-im, comm(F, G))          # L_F G = [F,G]/(i)

Random.seed!(20260825)
d = 3; M = 2
H = HS()
for m in 0:M; A = randn(ComplexF64,d,d); H[m] = m==0 ? (A+A')/2 : A; end
for m in 1:M; H[-m] = H[m]'; end

w = 1.0
tri = vanvleck(H, 3, w, d)                       # verified Deprit triangle, wd = 1

# Nikolaev eq (22) low orders, with H_2 = H_3 = ... = 0
H1 = H
nik0 = avg(Pk(H1,d), d)
nik1 = avg(scal(-1//2, Pk(Lk(H1, Sk(H1)), d)), d)
termA = scal(1//3, Pk(Lk(H1, Sk(Lk(H1, Sk(H1)))), d))
termB = scal(-1//6, Pk(Lk(H1, Sk(Sk(Lk(H1, Pk(H1,d))))), d))
nik2 = avg(plus(termA, termB), d)

println("Nikolaev explicit formula (eq 22 / §5.1) vs the verified Deprit triangle\n")
for (n, nk) in ((0,nik0),(1,nik1),(2,nik2))
    @printf("Heff^(%d):  ||triangle - Nikolaev|| = %.3e   (||triangle|| = %.3e)\n",
            n, opnorm(tri[n] - nk), opnorm(tri[n]))
end


L(G) = Lk(H1, G); S = Sk; P(G) = Pk(G, d)
add(xs...) = foldl(plus, xs)
# eps^4 line of H~ (H1-only terms), transcribed from the rendered page 12
nik3 = avg(add(
    scal( 1//6 , P(L(S(L(S(S(L(P(H1))))))))),
    scal(-1//4 , P(L(S(L(S(L(S(H1)))))))),
    scal( 1//12, P(L(S(S(L(S(L(P(H1))))))))),
    scal( 1//8 , P(L(S(S(L(P(L(S(H1))))))))),
    scal( 1//4 , P(L(P(L(S(S(L(S(H1))))))))),
    scal( 1//4 , P(L(P(L(S(L(S(S(H1))))))))),
    scal(-1//6 , P(L(P(L(S(S(S(L(P(H1)))))))))),
    scal(-1//4 , P(L(P(L(P(L(S(S(S(H1))))))))))), d)

tri4 = vanvleck(H, 4, w, d)

L(G) = Lk(H1, G); S = Sk; P(G) = Pk(G, d)
println("\nThe two identities the paper says it used (page 12 footnote):\n")
id1 = avg(P(L(S(H1))), d)          # P L_F S F  ==  0 ?
id2 = avg(P(L(P(H1))), d)          # P L_F P F  ==  0 ?
@printf("  ||P L_F S F|| = %.3e     <- this is exactly the eps^2 term of H~, i.e. Heff^(1)\n", opnorm(id1))
@printf("  ||P L_F P F|| = %.3e\n", opnorm(id2))
@printf("  ||Heff^(1)||  = %.3e     (from the verified triangle)\n", opnorm(vanvleck(H,2,w,d)[1]))
println()
println("  If 'P L_F S F == 0' held, the eps^2 term -(1/2) P L S H1 would vanish and Heff^(1)")
println("  would be zero. It is not (eq:Heff1 is nonzero). So that identity, as read, is not")
println("  what the paper means -> my reading of his notation is incomplete somewhere.")

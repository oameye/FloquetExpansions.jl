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

println("gauge check: van Vleck (Deprit triangle) vs canonical block-diagonalized Sambe\n")
for N in 1:4
    errs = Float64[]; ws = [40.0, 80.0, 160.0, 320.0]
    for w in ws
        hv = vanvleck(H, N, w, d)
        approx = sum(hv[n] for n in 0:N-1)
        exact  = sambe_heff(H, w, d, 14)
        push!(errs, opnorm(approx - exact))
    end
    sl = (log(errs[end]) - log(errs[1])) / (log(ws[end]) - log(ws[1]))
    @printf("N=%d  ||Heff^[N] - Heff_exact||: %s   slope = %+.3f  (expect %+.1f)\n",
            N, join([@sprintf("%.2e", e) for e in errs], "  "), sl, -N)
end

## The Deprit triangle. Notation and the derivation of the peeling recursion are in
## docs/dev/DERIVATION.md; the two families are named for their seed in docs/dev/DESIGN.md §4.

"""
    FloquetExpansion

Result of [`floquet_expansion`](@ref). Read it with [`effective_hamiltonian`](@ref) and
[`kick_operator`](@ref) rather than by field access.
"""
struct FloquetExpansion{G <: Gauge}
    H::PeriodicOperator
    K::Vector{PeriodicOperator}            # K[k] = K^(k), k = 1..N-1
    Kdot::Vector{PeriodicOperator}         # Kdot[k] = d/dt K^(k)
    dressedH::Vector{PeriodicOperator}     # ad_K^j on H,    flat over (n,j), UNWEIGHTED
    dressedKdot::Vector{PeriodicOperator}  # ad_K^j on Kdot, flat over (n,j), UNWEIGHTED
    Heff::Vector{SQA.QAdd}                 # Heff[n+1] = Heff^(n), n = 0..N-1
    order::Int
    wd::Symbolics.Num
end

# Bounds are known up front (n < N, j <= n), so the triangle is a flat vector rather than a
# Dict keyed by (n,j): statically sized and index-computable. See DESIGN.md §7 L1.
_tri(n::Int, j::Int) = (n * (n + 1)) ÷ 2 + j + 1

# Both families satisfy the SAME weight-free recursion, node(n,j) = sum_k [K^(k), node(n-k,j-1)],
# and differ only in seed; the Deprit weight is applied once, here, at assembly. Carrying the
# weight inside the recursion instead would multiply a rational per level, and SQA collapses a
# float-representable rational (1//2 -> 0.5 on the ComplexF64 tier) while keeping 1//3 exact, so
# the accumulated product silently drops to floats and leaves ~1e-16 residues in the result.
_weightH(j::Int) = im^j * (1 // factorial(j))
_weightKdot(j::Int) = im^j * (1 // factorial(j + 1))

Base.show(io::IO, ::MIME"text/plain", vv::FloquetExpansion{G}) where {G} =
    print(io, "FloquetExpansion{", nameof(G), "} of order ", vv.order)

# `wd` never appears in the bookkeeping: an order-n object carries wd^-n, so stripping it and
# tracking the power in the order index leaves the recursion wd-free. It is reattached on output.
_reattach(X::SQA.QAdd, wd, n::Int) = iszero(n) ? X : X * wd^(-n)
_reattach(X::PeriodicOperator, wd, n::Int) = iszero(n) ? X : wd^(-n) * X

# Operator-level zeros vanish for free, but coefficient-level cancellation on the symbolic tier
# does not happen without `simplify`, so mathematically-zero buckets survive as structurally
# nonzero and defeat pruning. Sweep once per order stage, on what propagates forward.
function _prune(X::PeriodicOperator)
    return PeriodicOperator(
        Dict{Int, SQA.QAdd}(l => SQA.simplify(Xl) for (l, Xl) in X.components),
    )
end

"""
    floquet_expansion(H::PeriodicOperator, wd, gauge::Gauge, order::Int) -> FloquetExpansion
    floquet_expansion(H::QAdd, w, t, gauge::Gauge, order::Int) -> FloquetExpansion

Expand a periodically driven Hamiltonian into a time-independent effective Hamiltonian and a
periodic kick operator, to the given `order` in `1/wd`.

The second form parses the time dependence first, so `H` may be given as an ordinary
time-dependent operator in the drive frequency `w` and time `t`.

Truncation follows the spec: `order = 1` is the rotating-wave approximation, and the error is
``\\mathcal{O}(w_d^{-\\text{order}})``. Note that the Liouvillian literature counts this shifted
by one.

The series is ASYMPTOTIC, not convergent. Beyond a problem-dependent optimal order, a higher
`order` makes the answer worse rather than better.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables w::Real t::Real g::Real;

julia> H = w * (a' * a) + g * cos(w * t) * (a + a');

julia> vv = floquet_expansion(H, w, t, VanVleck(), 1)
FloquetExpansion{VanVleck} of order 1

julia> effective_hamiltonian(vv)
w * a' * a
```

See also [`effective_hamiltonian`](@ref), [`kick_operator`](@ref), [`harmonics`](@ref).
"""
function floquet_expansion(H::PeriodicOperator, wd, gauge::Gauge, order::Int)
    order >= 1 || throw(ArgumentError("order must be >= 1, got $(order)"))
    LinearAlgebra.ishermitian(H) || throw(
        ArgumentError("the drive is not Hermitian: it must satisfy H_{-m} = H_m' (eq:fourierH)"),
    )

    nodes = (order * (order + 1)) ÷ 2
    dressedH = [zero(PeriodicOperator) for _ in 1:nodes]
    dressedKdot = [zero(PeriodicOperator) for _ in 1:nodes]
    K = PeriodicOperator[]
    Kdot = PeriodicOperator[]
    Heff = SQA.QAdd[]

    for n in 0:(order - 1)
        dressedH[_tri(n, 0)] = n == 0 ? H : zero(PeriodicOperator)

        for j in 1:n
            acc = zero(PeriodicOperator)
            for k in 1:(n - j + 1)
                acc = acc + SQA.commutator(K[k], dressedH[_tri(n - k, j - 1)])
            end
            dressedH[_tri(n, j)] = acc
        end

        for j in 1:n
            acc = zero(PeriodicOperator)
            for k in 1:(n - j + 1)
                # dressedKdot^(n)_[0] IS Kdot^(n+1) and is not computable at stage n. It is never
                # read either, since j >= 1 reaches only order n-k with k >= 1; read Kdot directly
                # so the unfillable slot is never touched.
                prev = j == 1 ? Kdot[n - k + 1] : dressedKdot[_tri(n - k, j - 1)]
                acc = acc + SQA.commutator(K[k], prev)
            end
            dressedKdot[_tri(n, j)] = acc
        end

        R = zero(PeriodicOperator)
        for j in 0:n
            R = R + _weightH(j) * dressedH[_tri(n, j)]
        end
        for j in 1:n
            R = R - _weightKdot(j) * dressedKdot[_tri(n, j)]
        end
        R = _prune(R)

        Heffn = SQA.simplify(time_average(R))
        push!(Heff, Heffn)

        # Stage order-1 would produce K^(order), which the approximant K^[N] = sum_{k<N} excludes
        # and no later stage reads. Skipping it saves the last antiderivative outright.
        if n < order - 1
            Knext = _prune(antiderivative(R - Heffn, gauge))
            push!(K, Knext)
            push!(Kdot, derivative(Knext))
        end
    end

    return FloquetExpansion{typeof(gauge)}(
        H, K, Kdot, dressedH, dressedKdot, Heff, order, Symbolics.Num(wd),
    )
end

function floquet_expansion(H::SQA.QAdd, w, t, gauge::Gauge, order::Int)
    return floquet_expansion(harmonics(H, w, t), w, gauge, order)
end

"""
    effective_hamiltonian(vv::FloquetExpansion) -> QAdd
    effective_hamiltonian(vv::FloquetExpansion, n::Int) -> QAdd

The time-independent effective Hamiltonian ``H_\\text{eff}^{[N]} = \\sum_{k<N} H_\\text{eff}^{(k)}``,
or with `n` the order-`n` contribution alone. The drive frequency is reattached, so the order-`n`
piece carries ``w_d^{-n}``.
"""
function effective_hamiltonian(vv::FloquetExpansion)
    out = zero(SQA.QAdd)
    for n in 0:(vv.order - 1)
        out = out + _reattach(vv.Heff[n + 1], vv.wd, n)
    end
    return SQA.simplify(out)
end

function effective_hamiltonian(vv::FloquetExpansion, n::Int)
    0 <= n < vv.order ||
        throw(ArgumentError("order $(n) is outside 0:$(vv.order - 1) for this expansion"))
    return SQA.simplify(_reattach(vv.Heff[n + 1], vv.wd, n))
end

"""
    kick_operator(vv::FloquetExpansion) -> PeriodicOperator
    kick_operator(vv::FloquetExpansion, t) -> QAdd

The periodic kick operator ``K^{[N]} = \\sum_{k<N} K^{(k)}`` generating the micromotion, as
harmonics or, given a time variable `t`, as the time-dependent operator ``K(t)``.

Under [`VanVleck`](@ref) this has vanishing time average, which is what makes
[`effective_hamiltonian`](@ref) independent of the initial phase of the drive.
"""
function kick_operator(vv::FloquetExpansion)
    out = zero(PeriodicOperator)
    for k in 1:length(vv.K)
        out = out + _reattach(vv.K[k], vv.wd, k)
    end
    return out
end

# Deliberately no `kick_operator(vv, n::Int)`: it would capture `kick_operator(vv, 0)`, which a
# caller means as t = 0, and silently return an order instead of a time.
kick_operator(vv::FloquetExpansion, t) = kick_operator(vv)(vv.wd, t)

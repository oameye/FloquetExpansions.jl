"""
    QuasienergyOperator(H::PeriodicOperator, nmax::Int)

The quasienergy operator ``Q = H_S - i\\partial_t`` in Sambe (extended) space, truncated to
harmonics `-nmax:nmax`:

```math
Q_{mn} = H_{m-n} - m\\,\\omega_d\\,\\delta_{mn}
```

Index it by harmonic, `Q[m, n]` with `m, n` in `-nmax:nmax`, not by array position.

The eigenvalues of `Q` are the quasienergies, defined modulo `H.wd`, and they are what
[`effective_hamiltonian`](@ref) approximates. This is an inspection and interop view: the
expansion itself reads harmonics directly and never materializes a Sambe matrix.

The sign of the diagonal follows the package's Fourier convention,
``H_S(t) = \\sum_m H_m e^{-i m \\omega_d t}``: acting on ``e^{-i n \\omega_d t}``,
``-i\\partial_t`` gives ``-n\\omega_d``. The opposite convention carries `+m*wd` and transposes
the off-diagonal to ``H_{n-m}``.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables w::Real t::Real;

julia> Q = QuasienergyOperator(
           harmonics(a' * a + a * expim(-w * t) + a' * expim(w * t), w, t), 1
       )
QuasienergyOperator over harmonics -1:1

julia> Q[1, 0]
a

julia> Q[1, 1]
-w + a' * a
```

See also [`PeriodicOperator`](@ref), [`floquet_expansion`](@ref).
"""
struct QuasienergyOperator
  blocks::Matrix{SQA.QAdd}
  nmax::Int
end

function QuasienergyOperator(H::PeriodicOperator, nmax::Int)
  nmax >= 0 || throw(ArgumentError("nmax must be >= 0, got $(nmax)"))
  n = 2nmax + 1
  blocks = Matrix{SQA.QAdd}(undef, n, n)
  for (i, m) in enumerate((-nmax):nmax), (j, k) in enumerate((-nmax):nmax)
    blocks[i, j] = m == k ? H[m - k] - (m * H.wd) * one(SQA.QAdd) : H[m - k]
  end
  return QuasienergyOperator(blocks, nmax)
end

function Base.getindex(Q::QuasienergyOperator, m::Int, n::Int)
  (abs(m) <= Q.nmax && abs(n) <= Q.nmax) || throw(BoundsError(Q, (m, n)))
  return Q.blocks[m + Q.nmax + 1, n + Q.nmax + 1]
end

Base.size(Q::QuasienergyOperator) = (2Q.nmax + 1, 2Q.nmax + 1)

"""
    harmonic_range(Q::QuasienergyOperator) -> UnitRange{Int}

The harmonics `Q` is truncated to, i.e. the valid index range for `Q[m, n]`.
"""
harmonic_range(Q::QuasienergyOperator) = (-Q.nmax):(Q.nmax)

function Base.show(io::IO, ::MIME"text/plain", Q::QuasienergyOperator)
  return print(io, "QuasienergyOperator over harmonics ", -Q.nmax, ":", Q.nmax)
end

Base.show(io::IO, Q::QuasienergyOperator) = show(io, MIME"text/plain"(), Q)

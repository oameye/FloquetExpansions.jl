"""
    QuasienergyOperator(G::PeriodicGenerator, nmax::Int)

Construct the quasienergy operator in Sambe (extended) space, truncated to harmonics
`-nmax:nmax`.
For a Hamiltonian generator, its blocks are

```math
Q_{mn} = H_{m-n} - m\\,\\omega_d\\,\\delta_{mn}
```

For a Liouvillian generator, the energy-like Floquet-Liouville operator has blocks

```math
Q_{mn} = i\\,\\mathcal{L}_{m-n} - m\\,\\omega_d\\,\\delta_{mn}.
```

Index it by harmonic, `Q[m, n]` with `m, n` in `-nmax:nmax`, rather than by array position.
The eigenvalues of `Q` are the quasienergies, defined modulo `G.wd`. For a dissipative
generator they are generally complex: the real part describes oscillation and the imaginary
part describes decay or growth. This type stores symbolic blocks only; numerical
vectorization is separate.

# Arguments

- `G`: Periodic Hamiltonian or Liouvillian generator.
- `nmax`: Maximum absolute harmonic label retained in the truncation.

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

julia> harmonic_range(Q)
-1:1
```

See also [`PeriodicGenerator`](@ref), [`floquet_expansion`](@ref).
"""
struct QuasienergyOperator{T}
  blocks::Matrix{T}
  nmax::Int
end

hamiltonian_quasienergy_component(component::SQA.QAdd) = component
liouvillian_quasienergy_component(component::Liouvillian) = im * component

function quasienergy_operator(
  G::PeriodicGenerator{T}, nmax::Int, identity::T, transform::F
) where {T,F}
  nmax >= 0 || throw(ArgumentError("nmax must be >= 0, got $(nmax)"))
  n = 2nmax + 1
  blocks = Matrix{T}(undef, n, n)
  for (i, m) in enumerate((-nmax):nmax), (j, k) in enumerate((-nmax):nmax)
    block = transform(G[m - k])
    blocks[i, j] = m == k ? block - (m * G.wd) * identity : block
  end
  return QuasienergyOperator{T}(blocks, nmax)
end

function QuasienergyOperator(G::PeriodicGenerator{SQA.QAdd}, nmax::Int)
  return quasienergy_operator(G, nmax, one(SQA.QAdd), hamiltonian_quasienergy_component)
end

function QuasienergyOperator(G::PeriodicGenerator{Liouvillian}, nmax::Int)
  return quasienergy_operator(G, nmax, one(Liouvillian), liouvillian_quasienergy_component)
end

function Base.getindex(Q::QuasienergyOperator, m::Int, n::Int)
  (abs(m) <= Q.nmax && abs(n) <= Q.nmax) || throw(BoundsError(Q, (m, n)))
  return Q.blocks[m + Q.nmax + 1, n + Q.nmax + 1]
end

Base.size(Q::QuasienergyOperator) = (2Q.nmax + 1, 2Q.nmax + 1)

"""
    harmonic_range(Q::QuasienergyOperator) -> UnitRange{Int}

Return the harmonic index range retained by `Q`, i.e. the valid index range for `Q[m, n]`.

See also [`QuasienergyOperator`](@ref).
"""
harmonic_range(Q::QuasienergyOperator) = (-Q.nmax):(Q.nmax)

function Base.show(io::IO, ::MIME"text/plain", Q::QuasienergyOperator)
  return print(io, "QuasienergyOperator over harmonics ", -Q.nmax, ":", Q.nmax)
end

Base.show(io::IO, Q::QuasienergyOperator) = show(io, MIME"text/plain"(), Q)

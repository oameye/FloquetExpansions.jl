const BOSON_SITE = 1
const PHASE_SITE = 2
const LEVEL_SITE = 3

struct AlgebraSite{R}
  kind::Int
  levels::Int
  hilbert::Vector{R}
end

boson_site(::Type{R}) where {R} = AlgebraSite{R}(BOSON_SITE, 0, R[])
phase_site(::Type{R}) where {R} = AlgebraSite{R}(PHASE_SITE, 0, R[])
function level_site(hilbert::Vector{R}) where {R}
  length(hilbert) >= 2 || throw(ArgumentError("a level site needs at least two levels"))
  all(>(0), hilbert) || throw(ArgumentError("a level site needs a positive Hilbert metric"))
  return AlgebraSite{R}(LEVEL_SITE, length(hilbert), hilbert)
end

const Monomial = Vector{Int}
const MonomialTerms{T} = Vector{Tuple{Monomial,T}}

struct OperatorAlgebra{T,R}
  sites::Vector{AlgebraSite{R}}
  products::Dict{Tuple{Monomial,Monomial},MonomialTerms{T}}
  adjoints::Dict{Monomial,MonomialTerms{T}}
end

function OperatorAlgebra{T}(sites::Vector{AlgebraSite{R}}) where {T,R}
  return OperatorAlgebra{T,R}(
    sites,
    Dict{Tuple{Monomial,Monomial},MonomialTerms{T}}(),
    Dict{Monomial,MonomialTerms{T}}(),
  )
end

identity_monomial(algebra::OperatorAlgebra) = zeros(Int, 2 * length(algebra.sites))

function site_identity_reduction(site::AlgebraSite, i::Int, j::Int, ::Type{T}) where {T}
  N = site.levels
  (i, j) == (N, N) || return Tuple{Tuple{Int,Int},T}[((i, j), one(T))]
  terms = Tuple{Tuple{Int,Int},T}[((0, 0), one(T))]
  for m in 1:(N - 1)
    push!(terms, ((m, m), -one(T)))
  end
  return terms
end

function site_product(
  site::AlgebraSite, left::Tuple{Int,Int}, right::Tuple{Int,Int}, ::Type{T}
) where {T}
  if site.kind == LEVEL_SITE
    left == (0, 0) && return Tuple{Tuple{Int,Int},T}[(right, one(T))]
    right == (0, 0) && return Tuple{Tuple{Int,Int},T}[(left, one(T))]
    left[2] == right[1] || return Tuple{Tuple{Int,Int},T}[]
    return site_identity_reduction(site, left[1], right[2], T)
  end
  k1, l1 = left
  k2, l2 = right
  terms = Tuple{Tuple{Int,Int},T}[]
  phase = site.kind == PHASE_SITE ? -im : one(T)
  for j in 0:min(l1, k2)
    c = T(binomial(big(l1), j) * binomial(big(k2), j) * factorial(big(j))) * phase^j
    push!(terms, ((k1 + k2 - j, l1 + l2 - j), c))
  end
  return terms
end

function site_adjoint(site::AlgebraSite, monomial::Tuple{Int,Int}, ::Type{T}) where {T}
  k, l = monomial
  if site.kind == LEVEL_SITE
    monomial == (0, 0) && return Tuple{Tuple{Int,Int},T}[((0, 0), one(T))]
    scale = T(site.hilbert[k] / site.hilbert[l])
    return [(m, scale * c) for (m, c) in site_identity_reduction(site, l, k, T)]
  elseif site.kind == BOSON_SITE
    return Tuple{Tuple{Int,Int},T}[((l, k), one(T))]
  end
  return site_product(site, (0, l), (k, 0), T)
end

function tensor_terms(algebra::OperatorAlgebra{T}, local_terms) where {T}
  result = MonomialTerms{T}([(Int[], one(T))])
  for terms in local_terms
    isempty(terms) && return MonomialTerms{T}()
    next = MonomialTerms{T}()
    for (monomial, c) in result, ((p, q), d) in terms
      push!(next, (vcat(monomial, [p, q]), c * d))
    end
    result = next
  end
  return result
end

function monomial_product(algebra::OperatorAlgebra{T}, A::Monomial, B::Monomial) where {T}
  return get!(algebra.products, (A, B)) do
    local_terms = [
      site_product(site, (A[2s - 1], A[2s]), (B[2s - 1], B[2s]), T) for
      (s, site) in pairs(algebra.sites)
    ]
    return tensor_terms(algebra, local_terms)
  end
end

function monomial_adjoint(algebra::OperatorAlgebra{T}, A::Monomial) where {T}
  return get!(algebra.adjoints, A) do
    local_terms = [
      site_adjoint(site, (A[2s - 1], A[2s]), T) for (s, site) in pairs(algebra.sites)
    ]
    return tensor_terms(algebra, local_terms)
  end
end

function monomial_degree(algebra::OperatorAlgebra, A::Monomial)
  degree = 0
  for (s, site) in pairs(algebra.sites)
    site.kind == LEVEL_SITE && continue
    degree += A[2s - 1] + A[2s]
  end
  return degree
end

struct AlgebraOperator{T,R}
  algebra::OperatorAlgebra{T,R}
  terms::Dict{Monomial,T}
end

struct AlgebraSuperoperator{T,R}
  algebra::OperatorAlgebra{T,R}
  terms::Dict{Tuple{Monomial,Monomial},T}
end

function accumulate!(terms::Dict{K,T}, key::K, value) where {K,T}
  iszero(value) && return terms
  updated = get(terms, key, zero(T)) + value
  iszero(updated) ? delete!(terms, key) : (terms[key] = updated)
  return terms
end

function algebra_operator(algebra::OperatorAlgebra{T,R}, pairs) where {T,R}
  terms = Dict{Monomial,T}()
  for (monomial, c) in pairs
    accumulate!(terms, Monomial(monomial), T(c))
  end
  return AlgebraOperator{T,R}(algebra, terms)
end

function Base.:+(A::AlgebraOperator{T,R}, B::AlgebraOperator{T,R}) where {T,R}
  terms = copy(A.terms)
  for (k, v) in B.terms
    accumulate!(terms, k, v)
  end
  return AlgebraOperator{T,R}(A.algebra, terms)
end

function Base.:*(c::Number, A::AlgebraOperator{T,R}) where {T,R}
  terms = Dict{Monomial,T}()
  for (k, v) in A.terms
    accumulate!(terms, k, T(c) * v)
  end
  return AlgebraOperator{T,R}(A.algebra, terms)
end

Base.:-(A::AlgebraOperator) = -1 * A
Base.:-(A::AlgebraOperator, B::AlgebraOperator) = A + (-B)
function Base.zero(A::AlgebraOperator{T,R}) where {T,R}
  return AlgebraOperator{T,R}(A.algebra, Dict{Monomial,T}())
end
Base.iszero(A::AlgebraOperator) = isempty(A.terms)
Base.:(==)(A::AlgebraOperator, B::AlgebraOperator) = A.terms == B.terms

function Base.:*(A::AlgebraOperator{T,R}, B::AlgebraOperator{T,R}) where {T,R}
  terms = Dict{Monomial,T}()
  for (a, ca) in A.terms, (b, cb) in B.terms
    for (m, c) in monomial_product(A.algebra, a, b)
      accumulate!(terms, m, ca * cb * c)
    end
  end
  return AlgebraOperator{T,R}(A.algebra, terms)
end

function Base.adjoint(A::AlgebraOperator{T,R}) where {T,R}
  terms = Dict{Monomial,T}()
  for (a, ca) in A.terms
    for (m, c) in monomial_adjoint(A.algebra, a)
      accumulate!(terms, m, conj(ca) * c)
    end
  end
  return AlgebraOperator{T,R}(A.algebra, terms)
end

function algebra_identity(algebra::OperatorAlgebra{T,R}) where {T,R}
  return AlgebraOperator{T,R}(algebra, Dict(identity_monomial(algebra) => one(T)))
end

function sandwich(A::AlgebraOperator{T,R}, B::AlgebraOperator{T,R}) where {T,R}
  terms = Dict{Tuple{Monomial,Monomial},T}()
  for (a, ca) in A.terms, (b, cb) in B.terms
    accumulate!(terms, (a, b), ca * cb)
  end
  return AlgebraSuperoperator{T,R}(A.algebra, terms)
end

left_action(A::AlgebraOperator) = sandwich(A, algebra_identity(A.algebra))
right_action(A::AlgebraOperator) = sandwich(algebra_identity(A.algebra), A)

function Base.:+(A::AlgebraSuperoperator{T,R}, B::AlgebraSuperoperator{T,R}) where {T,R}
  terms = copy(A.terms)
  for (k, v) in B.terms
    accumulate!(terms, k, v)
  end
  return AlgebraSuperoperator{T,R}(A.algebra, terms)
end

function Base.:*(c::Number, A::AlgebraSuperoperator{T,R}) where {T,R}
  terms = Dict{Tuple{Monomial,Monomial},T}()
  scale = T(c)
  iszero(scale) && return AlgebraSuperoperator{T,R}(A.algebra, terms)
  for (k, v) in A.terms
    terms[k] = scale * v
  end
  return AlgebraSuperoperator{T,R}(A.algebra, terms)
end

Base.:*(A::AlgebraSuperoperator, c::Number) = c * A
Base.:-(A::AlgebraSuperoperator) = -1 * A
Base.:-(A::AlgebraSuperoperator, B::AlgebraSuperoperator) = A + (-B)
Base.:(==)(A::AlgebraSuperoperator, B::AlgebraSuperoperator) = A.terms == B.terms
Base.iszero(A::AlgebraSuperoperator) = isempty(A.terms)
function Base.copy(A::AlgebraSuperoperator{T,R}) where {T,R}
  return AlgebraSuperoperator{T,R}(A.algebra, copy(A.terms))
end

function Base.zero(A::AlgebraSuperoperator{T,R}) where {T,R}
  return AlgebraSuperoperator{T,R}(A.algebra, Dict{Tuple{Monomial,Monomial},T}())
end

function Base.one(A::AlgebraSuperoperator{T,R}) where {T,R}
  identity = identity_monomial(A.algebra)
  return AlgebraSuperoperator{T,R}(A.algebra, Dict((identity, identity) => one(T)))
end

function Base.:*(A::AlgebraSuperoperator{T,R}, B::AlgebraSuperoperator{T,R}) where {T,R}
  algebra = A.algebra
  terms = Dict{Tuple{Monomial,Monomial},T}()
  for ((la, ra), ca) in A.terms, ((lb, rb), cb) in B.terms
    left = monomial_product(algebra, la, lb)
    right = monomial_product(algebra, rb, ra)
    for (l, cl) in left, (r, cr) in right
      accumulate!(terms, (l, r), ca * cb * cl * cr)
    end
  end
  return AlgebraSuperoperator{T,R}(algebra, terms)
end

struct ChargeBlock{F}
  entries::Vector{Tuple{Int,Int}}
  factor::F
end

struct ChargeGradedInverse{T,F} <: HomologicalInverse
  L0::Matrix{T}
  weights::Matrix{Int}
  regular::Vector{ChargeBlock{F}}
  singular::Vector{Vector{Int}}
  directions::Vector{Matrix{T}}
end

native_frame(representation::DenseLiouvilleRepresentation) = representation.basis
native_frame(representation::ExactLiouvilleRepresentation) = representation.frame

native_nonzero(x::Number, threshold::Real) = abs2(x) > threshold^2

function exact_nullspace(A::Matrix{R}) where {R}
  rows, columns = size(A)
  M = copy(A)
  pivots = Int[]
  row = 1
  for column in 1:columns
    row > rows && break
    pivot = exact_find_pivot(M, row, column)
    pivot == 0 && continue
    exact_swap_rows!(M, row, pivot)
    exact_scale_row!(M, row, inv(M[row, column]))
    for r in 1:rows
      factor = M[r, column]
      (r == row || iszero(factor)) && continue
      exact_row_axpy!(M, r, factor, row)
    end
    push!(pivots, column)
    row += 1
  end
  free = setdiff(1:columns, pivots)
  basis = Vector{R}[]
  for f in free
    v = exact_zeros(R, columns)
    v[f] = one(R)
    for (r, p) in pairs(pivots)
      v[p] = -M[r, f]
    end
    push!(basis, v)
  end
  return basis
end

function integer_vector(v::Vector{Rational{I}}) where {I<:Integer}
  scale = reduce(lcm, denominator.(v); init=one(I))
  scaled = v * scale
  divisor = reduce(gcd, numerator.(scaled); init=zero(I))
  iszero(divisor) && return zeros(Int, length(v))
  return Int[Int(numerator(x) ÷ divisor) for x in scaled]
end

function superoperator_indices(index::Int, d::Int)
  return (index - 1) % d + 1, (index - 1) ÷ d + 1
end

function charge_constraint!(rows, d, output, input)
  e, f = superoperator_indices(output, d)
  a, b = superoperator_indices(input, d)
  c = zeros(Int, d)
  c[e] += 1
  c[f] -= 1
  c[a] -= 1
  c[b] += 1
  iszero(c) || push!(rows, c)
  return rows
end

function charge_weights(L0::AbstractMatrix, operators, d::Int, threshold::Real)
  rows = Set{Vector{Int}}()
  for input in axes(L0, 2), output in axes(L0, 1)
    native_nonzero(L0[output, input], threshold) &&
      charge_constraint!(rows, d, output, input)
  end
  for X in operators
    support = [(a, b) for b in 1:d, a in 1:d if native_nonzero(X[a, b], threshold)]
    for k in 2:length(support)
      (a, b), (c, e) = support[1], support[k]
      constraint = zeros(Int, d)
      constraint[a] += 1
      constraint[b] -= 1
      constraint[c] -= 1
      constraint[e] += 1
      iszero(constraint) || push!(rows, constraint)
    end
  end
  gauge = zeros(Int, d)
  gauge[1] = 1
  push!(rows, gauge)
  A = Rational{BigInt}[row[j] for row in collect(rows), j in 1:d]
  kernel = exact_nullspace(A)
  weights = zeros(Int, length(kernel), d)
  for (r, v) in pairs(kernel)
    weights[r, :] = integer_vector(v)
  end
  return weights
end

function superoperator_label(weights::Matrix{Int}, d::Int, output::Int, input::Int)
  e, f = superoperator_indices(output, d)
  a, b = superoperator_indices(input, d)
  return weights[:, e] - weights[:, f] - weights[:, a] + weights[:, b]
end

function charge_sectors(weights::Matrix{Int}, d::Int)
  sectors = Dict{Vector{Int},Vector{Tuple{Int,Int}}}()
  for input in 1:(d ^ 2), output in 1:(d ^ 2)
    label = superoperator_label(weights, d, output, input)
    push!(get!(sectors, label, Tuple{Int,Int}[]), (output, input))
  end
  return sectors
end

function adjoint_block(L0::AbstractMatrix{T}, entries::Vector{Tuple{Int,Int}}) where {T}
  m = length(entries)
  A = exact_zeros(T, m, m)
  position = Dict(entry => k for (k, entry) in pairs(entries))
  for (column, (output, input)) in pairs(entries)
    for row in axes(L0, 1)
      k = get(position, (row, input), 0)
      k == 0 || (A[k, column] += L0[row, output])
    end
    for col in axes(L0, 2)
      k = get(position, (output, col), 0)
      k == 0 || (A[k, column] -= L0[input, col])
    end
  end
  return A
end

function block_factorization(A::Matrix{ComplexF64}, tol::Real)
  values = LinearAlgebra.svdvals(A)
  minimum(values) > tol * max(1.0, maximum(values)) || return nothing
  return LinearAlgebra.lu(A)
end

function block_factorization(A::Matrix{T}, ::Real) where {T}
  length(exact_row_basis(A)) == size(A, 1) || return nothing
  return exact_solve(A, Matrix{T}(LinearAlgebra.I, size(A)...))
end

block_solve(factor::LinearAlgebra.LU, b) = factor \ b
block_solve(inverse::Matrix, b) = inverse * b

function independent_directions(projected::Vector{Matrix{ComplexF64}}, tol::Real)
  isempty(projected) && return projected
  M = reduce(hcat, [vcat(real(vec(G)), imag(vec(G))) for G in projected])
  F = LinearAlgebra.svd(M)
  rank = count(>(tol * max(1.0, F.S[1])), F.S)
  return Matrix{ComplexF64}[
    sum(F.V[j, i] * projected[j] for j in eachindex(projected)) / F.S[i] for i in 1:rank
  ]
end

function independent_directions(projected::Vector{Matrix{T}}, ::Real) where {T}
  nonzero = filter(!iszero, projected)
  isempty(nonzero) && return nonzero
  coordinates = reduce(
    vcat, [transpose(vcat(real.(vec(G)), imag.(vec(G)))) for G in nonzero]
  )
  return nonzero[exact_row_basis(coordinates)]
end

function ChargeGradedInverse(
  representation::NativeRepresentation, L0::Matrix{T}, leading::AbstractMatrix, tol::Real
) where {T}
  d = round(Int, sqrt(size(L0, 1)))
  d^2 == size(L0, 1) || throw(DimensionMismatch("L0 must act on a d^2-dimensional space"))
  frame = native_frame(representation)
  operators = [
    sum(leading[k, j] * frame[k] for k in eachindex(frame)) for j in axes(leading, 2)
  ]
  threshold = T === ComplexF64 ? tol * max(1.0, LinearAlgebra.norm(L0)) : 0.0
  weights = charge_weights(L0, operators, d, threshold)
  sectors = charge_sectors(weights, d)

  first_factor = block_factorization(exact_zeros(T, 1, 1) + LinearAlgebra.I, tol)
  regular = ChargeBlock{typeof(first_factor)}[]
  singular = Vector{Int}[]
  for (label, entries) in sectors
    iszero(label) && (push!(singular, label); continue)
    factor = block_factorization(adjoint_block(L0, entries), tol)
    if factor === nothing
      push!(singular, label)
    else
      push!(regular, ChargeBlock(entries, factor))
    end
  end

  keep = falses(d^2, d^2)
  for label in singular, (output, input) in sectors[label]
    keep[output, input] = true
  end
  projected = [G .* keep for G in native_gauge_directions(representation)]
  directions = independent_directions(projected, tol)
  return ChargeGradedInverse(Matrix{T}(L0), weights, regular, singular, directions)
end

function regular_gauge(
  inverse::ChargeGradedInverse{T},
  representation::NativeRepresentation,
  L0,
  residual,
  known,
  active,
  tol,
) where {T}
  native_matches(representation, L0, inverse.L0, tol) || throw(
    ArgumentError("charge-graded inverse was built for a different averaged generator")
  )
  S = zero(inverse.L0)
  isempty(inverse.regular) && return S
  target = native_dark_target(representation, residual, known, active, tol)
  d = size(target, 1) == 0 ? 0 : round(Int, sqrt(size(L0, 1)))
  Y = native_gksl(representation, exact_zeros(T, d, d), target)
  for block in inverse.regular
    x = block_solve(block.factor, T[-Y[output, input] for (output, input) in block.entries])
    for (k, (output, input)) in pairs(block.entries)
      S[output, input] = x[k]
    end
  end
  return S
end

function singular_gauge_directions(inverse::ChargeGradedInverse, ::NativeRepresentation)
  return inverse.directions
end

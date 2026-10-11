const ExactReal = Rational{BigInt}
const ExactPoint = Vector{ExactReal}
const SparsePolynomial = Vector{Tuple{Vector{Int},ExactReal}}

const FIT_MARGIN = 2
const CLASS_POINTS = 6
const MAXIMUM_LINE_NODES = 24
const MAXIMUM_UNKNOWNS = 250
const MAXIMUM_FIT_ROUNDS = 6
const VERIFICATION_SPREAD = 3
const PRIME_RATIOS = (101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163)

function reconstruction_error(message::String)
  return ArgumentError(
    "exact polynomial reconstruction of the symbolic coefficients failed: $message. " *
    "Substitute exact values for some parameters instead.",
  )
end

function monomial_value(exponents::Vector{Int}, x::ExactPoint)
  value = one(ExactReal)
  for (xi, e) in zip(x, exponents)
    iszero(e) || (value *= xi^e)
  end
  return value
end

function polynomial_value(p::SparsePolynomial, x::ExactPoint)
  value = zero(ExactReal)
  for (exponents, c) in p
    value += c * monomial_value(exponents, x)
  end
  return value
end

function monomial_support(box::Vector{Int}, low::Int, high::Int)
  support = Vector{Int}[]
  for index in CartesianIndices(Tuple(0:d for d in box))
    exponents = collect(Tuple(index))
    low <= sum(exponents) <= high && push!(support, exponents)
  end
  return support
end

function eliminate!(A::Matrix{ExactReal}, columns::Int)
  pivots = Int[]
  for column in 1:columns
    row = length(pivots) + 1
    row > size(A, 1) && break
    pivot = findfirst(i -> !iszero(A[i, column]), row:size(A, 1))
    pivot === nothing && continue
    pivot += row - 1
    for j in axes(A, 2)
      A[row, j], A[pivot, j] = A[pivot, j], A[row, j]
    end
    scale = A[row, column]
    for j in axes(A, 2)
      A[row, j] /= scale
    end
    for i in axes(A, 1)
      (i == row || iszero(A[i, column])) && continue
      factor = A[i, column]
      for j in axes(A, 2)
        A[i, j] -= factor * A[row, j]
      end
    end
    push!(pivots, column)
  end
  return pivots
end

function line_polynomial(xs::Vector{ExactReal}, ys::Vector{ExactReal})
  n = length(xs)
  c = copy(ys)
  for j in 2:n, i in n:-1:j
    c[i] = (c[i] - c[i - 1]) / (xs[i] - xs[i - j + 1])
  end
  p = ExactReal[c[n]]
  for i in (n - 1):-1:1
    q = zeros(ExactReal, length(p) + 1)
    q[2:end] .+= p
    q[1:(end - 1)] .-= xs[i] .* p
    q[1] += c[i]
    p = q
  end
  degree = something(findlast(!iszero, p), 1)
  return degree + FIT_MARGIN < n ? p[1:degree] : nothing
end

struct Sampler{K,F}
  evaluate::F
  signs::Vector{Int}
  points::Vector{ExactPoint}
  samples::Vector{Dict{K,ExactReal}}
end

function Sampler{K}(evaluate::F, signs::Vector{Int}) where {K,F}
  return Sampler{K,F}(evaluate, signs, ExactPoint[], Dict{K,ExactReal}[])
end

function sample!(sampler::Sampler, magnitudes::ExactPoint)
  point = sampler.signs .* magnitudes
  index = findfirst(isequal(point), sampler.points)
  index === nothing || return index
  push!(sampler.samples, sampler.evaluate(point))
  push!(sampler.points, point)
  return length(sampler.points)
end

function pseudo_random(state::Int)
  return mod(state * 1103515245 + 12345, 2^31)
end

function random_magnitudes(m::Int, seed::Int, spread::Int)
  state = seed
  magnitudes = ExactReal[]
  for _ in 1:m
    draws = Int[]
    for _ in 1:3
      state = pseudo_random(state)
      push!(draws, state >> 8)
    end
    scale = ExactReal(3)^(draws[3] % (2spread + 1) - spread)
    push!(magnitudes, (1 + draws[1] % 9) // (2 + draws[2] % 8) * scale)
  end
  return magnitudes
end

function reference_magnitudes(m::Int)
  n = length(PRIME_RATIOS)
  return ExactReal[PRIME_RATIOS[mod1(2i, n)] // PRIME_RATIOS[mod1(3i + 1, n)] for i in 1:m]
end

function value_classes(sampler::Sampler{K}, keys::Vector{K}) where {K}
  classes = Dict{Vector{ExactReal},Vector{K}}()
  for key in keys
    column = ExactReal[get(sample, key, zero(ExactReal)) for sample in sampler.samples]
    push!(get!(() -> K[], classes, column), key)
  end
  return [sort!(members; by=string) for members in values(classes)]
end

function scan_line(sampler::Sampler, keys::Vector, line)
  xs, indices = ExactReal[], Int[]
  found = Dict{Any,Vector{ExactReal}}()
  for j in 1:MAXIMUM_LINE_NODES
    push!(xs, j // (j + 1))
    push!(indices, sample!(sampler, line(xs[end])))
    for key in keys
      haskey(found, key) && continue
      y = ExactReal[get(sampler.samples[i], key, zero(ExactReal)) for i in indices]
      p = line_polynomial(xs, y)
      p === nothing || (found[key] = p)
    end
    length(found) == length(keys) && return found
  end
  return throw(
    reconstruction_error(
      "a coefficient is not a polynomial of degree below $(MAXIMUM_LINE_NODES - FIT_MARGIN - 1) " *
      "in the parameters",
    ),
  )
end

function support_models(sampler::Sampler, keys::Vector, base::ExactPoint)
  m = length(base)
  boxes = [scan_line(sampler, keys, t -> setindex!(copy(base), t, i)) for i in 1:m]
  rays = scan_line(sampler, keys, t -> t .* base)
  return Dict(
    key => monomial_support(
      [length(boxes[i][key]) - 1 for i in 1:m],
      something(findfirst(!iszero, rays[key]), 1) - 1,
      length(rays[key]) - 1,
    ) for key in keys
  )
end

function fit_group(sampler::Sampler, rows::Vector{Int}, keys::Vector, support)
  u = length(support)
  rows = rows[(end - u - FIT_MARGIN + 1):end]
  A = hcat(
    ExactReal[monomial_value(e, sampler.points[i]) for i in rows, e in support],
    ExactReal[get(sampler.samples[i], key, zero(ExactReal)) for i in rows, key in keys],
  )
  length(eliminate!(A, u)) == u || return nothing
  all(iszero, A[(u + 1):end, (u + 1):end]) || return nothing
  return Dict(
    key => Tuple{Vector{Int},ExactReal}[
      (e, c) for (e, c) in zip(support, A[1:u, u + k]) if !iszero(c)
    ] for (k, key) in pairs(keys)
  )
end

function fit_models(sampler::Sampler, models::Dict)
  groups = Dict{Vector{Vector{Int}},Vector{Any}}()
  for (key, support) in models
    push!(get!(() -> Any[], groups, support), key)
  end
  needed = maximum(length, keys(groups))
  needed <= MAXIMUM_UNKNOWNS || throw(
    reconstruction_error(
      "a coefficient has more than $MAXIMUM_UNKNOWNS candidate monomials"
    ),
  )
  fits = Dict{Any,SparsePolynomial}()
  rows = Int[]
  seed = 0
  for round in 1:MAXIMUM_FIT_ROUNDS
    while length(rows) < needed + FIT_MARGIN + 2round
      seed += 1
      push!(rows, sample!(sampler, random_magnitudes(length(sampler.signs), seed, 1)))
    end
    for (support, members) in groups
      all(key -> haskey(fits, key), members) && continue
      merge!(fits, something(fit_group(sampler, rows, members, support), Dict()))
    end
    length(fits) == length(models) && return fits
  end
  return throw(
    reconstruction_error("the exact samples fit no polynomial of the detected degrees")
  )
end

function verify_fits(sampler::Sampler, functions::Dict, count::Int)
  verified = ExactPoint[]
  for v in 1:count
    magnitudes = random_magnitudes(length(sampler.signs), 5000 + v, VERIFICATION_SPREAD)
    index = sample!(sampler, magnitudes)
    point, sample = sampler.points[index], sampler.samples[index]
    for (key, f) in functions
      polynomial_value(f, point) == get(sample, key, zero(ExactReal)) || throw(
        reconstruction_error(
          "a reconstructed coefficient disagrees with the exact solver at the point " *
          "$(join(string.(point), ", "))",
        ),
      )
    end
    push!(verified, point)
  end
  return verified
end

function reconstruct_polynomials(
  evaluate, keys::Vector{K}, signs::Vector{Int}, verification::Int
) where {K}
  m = length(signs)
  sampler = Sampler{K}(evaluate, signs)
  base = reference_magnitudes(m)
  sample!(sampler, base)
  for seed in 1:CLASS_POINTS
    sample!(sampler, random_magnitudes(m, 1000 + seed, 1))
  end
  classes = value_classes(sampler, keys)
  fits = fit_models(sampler, support_models(sampler, first.(classes), base))
  functions = Dict(key => fits[first(class)] for class in classes for key in class)
  return functions, verify_fits(sampler, functions, verification)
end

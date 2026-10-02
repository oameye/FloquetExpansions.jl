const ExactField = Complex{Rational{BigInt}}
const FieldPolynomial = Vector{ExactField}

function poly_degree(p::FieldPolynomial)
  last = findlast(!iszero, p)
  return last === nothing ? -1 : last - 1
end

poly_trim(p::FieldPolynomial) = p[1:(poly_degree(p) + 1)]

function poly_valuation(p::FieldPolynomial)
  first = findfirst(!iszero, p)
  return first === nothing ? -1 : first - 1
end

function poly_sub(p::FieldPolynomial, q::FieldPolynomial)
  result = zeros(ExactField, max(length(p), length(q)))
  result[eachindex(p)] .= p
  for (k, c) in pairs(q)
    result[k] -= c
  end
  return poly_trim(result)
end

function poly_mul(p::FieldPolynomial, q::FieldPolynomial)
  (isempty(p) || isempty(q)) && return ExactField[]
  result = zeros(ExactField, length(p) + length(q) - 1)
  for (j, b) in pairs(q), (i, a) in pairs(p)
    result[i + j - 1] += a * b
  end
  return poly_trim(result)
end

function poly_eval(p::FieldPolynomial, x::ExactField)
  value = zero(ExactField)
  for k in reverse(eachindex(p))
    value = value * x + p[k]
  end
  return value
end

function poly_divrem(a::FieldPolynomial, b::FieldPolynomial)
  db = poly_degree(b)
  r = copy(a)
  q = zeros(ExactField, max(poly_degree(a) - db + 1, 0))
  for k in (poly_degree(a) - db):-1:0
    c = r[k + db + 1] / b[db + 1]
    q[k + 1] = c
    iszero(c) && continue
    for j in 0:db
      r[k + j + 1] -= c * b[j + 1]
    end
  end
  return poly_trim(q), poly_trim(r)
end

function node_polynomial(xs::Vector{ExactField})
  p = ExactField[one(ExactField)]
  for x in xs
    p = poly_mul(p, ExactField[-x, one(ExactField)])
  end
  return p
end

function interpolating_polynomial(xs::Vector{ExactField}, ys::Vector{ExactField})
  n = length(xs)
  c = copy(ys)
  for j in 2:n, i in n:-1:j
    c[i] = (c[i] - c[i - 1]) / (xs[i] - xs[i - j + 1])
  end
  p = ExactField[c[n]]
  for i in (n - 1):-1:1
    p = poly_sub(poly_mul(p, ExactField[-xs[i], one(ExactField)]), ExactField[-c[i]])
  end
  return poly_trim(p)
end

function rational_reconstruction(
  xs::Vector{ExactField}, ys::Vector{ExactField}, margin::Int
)
  r0, r1 = node_polynomial(xs), interpolating_polynomial(xs, ys)
  t0, t1 = ExactField[], ExactField[one(ExactField)]
  best = nothing
  gap = margin
  while poly_degree(r1) >= 0
    q, r = poly_divrem(r0, r1)
    poly_degree(q) > gap && ((gap, best) = (poly_degree(q), (r1, t1)))
    r0, r1 = r1, r
    t0, t1 = t1, poly_sub(t0, poly_mul(q, t1))
  end
  best === nothing && return nothing
  numerator, denominator = best
  all(x -> !iszero(poly_eval(denominator, x)), xs) || return nothing
  return numerator, denominator
end

struct RationalFunction
  numerator::Vector{Tuple{Vector{Int},ExactField}}
  denominator::Vector{Tuple{Vector{Int},ExactField}}
end

function Base.:(==)(f::RationalFunction, g::RationalFunction)
  return f.numerator == g.numerator && f.denominator == g.denominator
end

Base.hash(f::RationalFunction, h::UInt) = hash(f.denominator, hash(f.numerator, h))

function monomial_value(exponents::Vector{Int}, point::Vector{ExactField})
  value = one(ExactField)
  for (x, e) in zip(point, exponents)
    value *= x^e
  end
  return value
end

function sparse_value(
  terms::Vector{Tuple{Vector{Int},ExactField}}, point::Vector{ExactField}
)
  value = zero(ExactField)
  for (exponents, c) in terms
    value += c * monomial_value(exponents, point)
  end
  return value
end

function rational_value(f::RationalFunction, point::Vector{ExactField})
  return sparse_value(f.numerator, point) / sparse_value(f.denominator, point)
end

function monomial_support(box::Vector{Int}, low::Int, high::Int)
  support = Vector{Int}[]
  for index in CartesianIndices(Tuple(0:d for d in box))
    exponents = collect(Tuple(index))
    low <= sum(exponents; init=0) <= high && push!(support, exponents)
  end
  return support
end

function evaluation_matrix(support::Vector{Vector{Int}}, points::Vector{Vector{ExactField}})
  return ExactField[monomial_value(e, x) for x in points, e in support]
end

function eliminate_column!(A::Matrix{ExactField}, row::Int, column::Int)
  pivot = findfirst(i -> !iszero(A[i, column]), row:size(A, 1))
  pivot === nothing && return false
  pivot += row - 1
  A[row, :], A[pivot, :] = A[pivot, :], A[row, :]
  A[row, :] ./= A[row, column]
  for i in axes(A, 1)
    (i == row || iszero(A[i, column])) && continue
    A[i, :] .-= A[i, column] .* A[row, :]
  end
  return true
end

function row_reduce!(A::Matrix{ExactField}, columns::Int)
  pivots = Int[]
  for column in 1:columns
    length(pivots) == size(A, 1) && break
    eliminate_column!(A, length(pivots) + 1, column) && push!(pivots, column)
  end
  return pivots
end

struct ExactSampler{K,F}
  evaluate::F
  keys::Vector{K}
  signs::Vector{Int}
  points::Vector{Vector{Rational{BigInt}}}
  samples::Vector{Dict{K,ExactField}}
end

function ExactSampler(evaluate, keys::Vector{K}, signs::Vector{Int}) where {K}
  return ExactSampler(
    evaluate, keys, signs, Vector{Rational{BigInt}}[], Dict{K,ExactField}[]
  )
end

function sample!(sampler::ExactSampler, magnitudes::Vector{Rational{BigInt}})
  values = sampler.signs .* magnitudes
  index = findfirst(isequal(values), sampler.points)
  index === nothing || return sampler.samples[index]
  sample = sampler.evaluate(values)
  push!(sampler.points, values)
  push!(sampler.samples, sample)
  return sample
end

function small_rational(seed::Int)
  first = mod(seed * 1103515245 + 12345, 2^31)
  second = mod(first * 1103515245 + 12345, 2^31)
  return Rational{BigInt}(1 + mod(first >> 8, 9), 2 + mod(second >> 8, 8))
end

random_values(m::Int, k::Int) = Rational{BigInt}[small_rational(97k + i) for i in 1:m]

function coprime_numerators(denominator::Int)
  return [p for p in 1:(2denominator - 1) if gcd(p, denominator) == 1]
end

function line_node(j::Int)
  denominator = 2
  while true
    numerators = coprime_numerators(denominator)
    j <= length(numerators) && return Rational{BigInt}(numerators[j], denominator)
    j -= length(numerators)
    denominator += 1
  end
end

const RECONSTRUCTION_MARGIN = 2
const MAXIMUM_LINE_SAMPLES = 32
const MAXIMUM_UNKNOWNS = 400
const MAXIMUM_FIT_ROUNDS = 8
const VERIFICATION_POINTS = 2

const LineFit = Tuple{FieldPolynomial,FieldPolynomial}

function reconstruction_error(message::String)
  return ArgumentError(
    "exact rational reconstruction of the symbolic coefficients failed: $message. " *
    "Substitute exact values for some parameters instead.",
  )
end

function record_line_fits!(fits::Dict, pending::Set, xs, ys::Dict)
  for key in collect(pending)
    fit = rational_reconstruction(xs, ys[key], RECONSTRUCTION_MARGIN)
    fit === nothing && continue
    fits[key] = fit
    delete!(pending, key)
  end
  return fits
end

function line_fits(sampler::ExactSampler{K}, line) where {K}
  pending = Set{K}(sampler.keys)
  xs = ExactField[]
  ys = Dict{K,Vector{ExactField}}(key => ExactField[] for key in sampler.keys)
  fits = Dict{K,LineFit}()
  for j in 1:MAXIMUM_LINE_SAMPLES
    x = line_node(j)
    sample = sample!(sampler, line(x))
    push!(xs, ExactField(x))
    foreach(key -> push!(ys[key], sample[key]), sampler.keys)
    record_line_fits!(fits, pending, xs, ys)
    isempty(pending) && return fits
  end
  return throw(
    reconstruction_error("a coefficient needs more than $MAXIMUM_LINE_SAMPLES line samples")
  )
end

function coordinate_line(base::Vector{Rational{BigInt}}, i::Int)
  return function (x)
    point = copy(base)
    point[i] = x
    return point
  end
end

struct OutputModel
  numerator::Vector{Vector{Int}}
  denominator::Vector{Vector{Int}}
end

function is_polynomial_model(model::OutputModel)
  return length(model.denominator) == 1 && all(iszero, model.denominator[1])
end

function variable_degrees(lines::Vector, keys_::Vector)
  m = length(lines)
  return Dict(
    key => (
      [poly_degree(lines[i][key][1]) for i in 1:m],
      [poly_degree(lines[i][key][2]) for i in 1:m],
    ) for key in keys_
  )
end

function ray_fits(sampler::ExactSampler, lines::Vector, base::Vector{Rational{BigInt}})
  length(base) == 1 && return only(lines)
  return line_fits(sampler, t -> t .* base)
end

function output_models(sampler::ExactSampler{K}, base::Vector{Rational{BigInt}}) where {K}
  lines = [line_fits(sampler, coordinate_line(base, i)) for i in eachindex(base)]
  degrees = variable_degrees(lines, sampler.keys)
  ray = ray_fits(sampler, lines, base)
  rational = Set{K}(
    key for (key, (_, d)) in degrees if any(>(0), d) || poly_degree(ray[key][2]) > 0
  )
  direction = random_values(length(base), 0)
  shifted = if isempty(rational) || length(base) == 1
    ray
  else
    line_fits(sampler, t -> base .+ t .* direction)
  end
  return Dict{K,OutputModel}(
    key => output_model(degrees[key], ray[key], shifted[key], key in rational) for
    key in sampler.keys
  )
end

function output_model(degrees, ray::LineFit, shifted::LineFit, rational::Bool)
  numerator, denominator = degrees
  rational || return OutputModel(
    monomial_support(numerator, poly_valuation(ray[1]), poly_degree(ray[1])),
    [zeros(Int, length(numerator))],
  )
  return OutputModel(
    monomial_support(numerator, 0, poly_degree(shifted[1])),
    monomial_support(denominator, 0, poly_degree(shifted[2])),
  )
end

sample_points(sampler::ExactSampler) = [ExactField.(values) for values in sampler.points]

function inconsistent_samples_error()
  return reconstruction_error(
    "the exact samples are not consistent with a rational function of the detected degrees"
  )
end

function sparse_terms(support::Vector{Vector{Int}}, coefficients::AbstractVector)
  return [(e, ExactField(c)) for (e, c) in zip(support, coefficients) if !iszero(c)]
end

function polynomial_fit(support::Vector{Vector{Int}}, keys_::Vector, sampler::ExactSampler)
  u = length(support)
  values = ExactField[sample[key] for sample in sampler.samples, key in keys_]
  A = hcat(evaluation_matrix(support, sample_points(sampler)), values)
  pivots = row_reduce!(A, u)
  length(pivots) < u && return u - length(pivots), RationalFunction[]
  all(iszero, A[(u + 1):end, (u + 1):end]) || throw(inconsistent_samples_error())
  one_term = [(zeros(Int, length(sampler.signs)), one(ExactField))]
  functions = [
    RationalFunction(sparse_terms(support, A[1:u, u + k]), one_term) for
    k in eachindex(keys_)
  ]
  return 0, functions
end

function nullspace_vector(A::Matrix{ExactField}, pivots::Vector{Int})
  free = only(setdiff(axes(A, 2), pivots))
  v = zeros(ExactField, size(A, 2))
  v[free] = one(ExactField)
  for (row, column) in pairs(pivots)
    v[column] = -A[row, free]
  end
  return v
end

function rational_fit(model::OutputModel, key, sampler::ExactSampler)
  points = sample_points(sampler)
  y = ExactField[sample[key] for sample in sampler.samples]
  A = hcat(
    evaluation_matrix(model.numerator, points),
    -y .* evaluation_matrix(model.denominator, points),
  )
  pivots = row_reduce!(A, size(A, 2))
  deficit = size(A, 2) - 1 - length(pivots)
  deficit > 0 && return deficit, RationalFunction[]
  deficit < 0 && throw(inconsistent_samples_error())
  v = nullspace_vector(A, pivots)
  u = length(model.numerator)
  lead = findfirst(!iszero, v[(u + 1):end])
  lead === nothing && throw(inconsistent_samples_error())
  v ./= v[u + lead]
  f = RationalFunction(
    sparse_terms(model.numerator, v[1:u]), sparse_terms(model.denominator, v[(u + 1):end])
  )
  return 0, [f]
end

function model_groups(models::Dict{K,OutputModel}) where {K}
  groups = Dict{Vector{Vector{Int}},Vector{K}}()
  for (key, model) in models
    is_polynomial_model(model) || continue
    push!(get!(() -> K[], groups, model.numerator), key)
  end
  return groups
end

function fit_round(sampler::ExactSampler{K}, models::Dict{K,OutputModel}) where {K}
  fitted = Dict{K,RationalFunction}()
  deficit = 0
  for (support, keys_) in model_groups(models)
    d, functions = polynomial_fit(support, keys_, sampler)
    deficit = max(deficit, d)
    d == 0 && merge!(fitted, Dict(zip(keys_, functions)))
  end
  for (key, model) in models
    is_polynomial_model(model) && continue
    d, functions = rational_fit(model, key, sampler)
    deficit = max(deficit, d)
    d == 0 && (fitted[key] = only(functions))
  end
  return deficit, fitted
end

function check_unknowns(models::Dict)
  for model in values(models)
    length(model.numerator) + length(model.denominator) <= MAXIMUM_UNKNOWNS || throw(
      reconstruction_error(
        "a coefficient has more than $MAXIMUM_UNKNOWNS candidate monomials"
      ),
    )
  end
  return models
end

function fit_outputs(sampler::ExactSampler, models::Dict)
  check_unknowns(models)
  m = length(sampler.signs)
  seed = 1
  for _ in 1:MAXIMUM_FIT_ROUNDS
    deficit, fitted = fit_round(sampler, models)
    deficit == 0 && return fitted
    for _ in 1:deficit
      sample!(sampler, random_values(m, seed))
      seed += 1
    end
  end
  return throw(inconsistent_samples_error())
end

function verify_reconstruction(
  fitted::Dict{K,RationalFunction}, values::Vector{Rational{BigInt}}, sample::Dict{K}
) where {K}
  point = ExactField.(values)
  for (key, value) in sample
    rational_value(fitted[key], point) == value || throw(
      reconstruction_error(
        "the reconstruction does not match the exact samples at $(join(string.(values), ", "))",
      ),
    )
  end
  return fitted
end

function verify_fresh_points(sampler::ExactSampler, fitted::Dict)
  m = length(sampler.signs)
  verified = 0
  seed = 1000
  while verified < VERIFICATION_POINTS
    values = sampler.signs .* random_values(m, seed)
    seed += 1
    any(isequal(values), sampler.points) && continue
    verify_reconstruction(fitted, values, sampler.evaluate(values))
    verified += 1
  end
  return fitted
end

function reconstruct_rational_functions(sampler::ExactSampler)
  models = output_models(sampler, random_values(length(sampler.signs), -1))
  return verify_fresh_points(sampler, fit_outputs(sampler, models))
end

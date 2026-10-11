const TimeMonomial = Tuple{Vector{Int},Vector{Int}}
const KernelKey = Tuple{Int,TimeMonomial}

struct TimeKernel{X}
  terms::Dict{KernelKey,X}
end

function kernel_identity(identity::X, q::Int) where {X}
  return TimeKernel{X}(Dict{KernelKey,X}((0, (zeros(Int, q), zeros(Int, q))) => identity))
end

function shifted_monomial(
  (phases, powers)::TimeMonomial, variable::Int, phase::Int, power::Int
)
  shifted_phases = copy(phases)
  shifted_powers = copy(powers)
  shifted_phases[variable] += phase
  shifted_powers[variable] += power
  return (shifted_phases, shifted_powers)
end

function difference_expansion(p::Int, upper::Int, lower::Int, q::Int)
  upper == q + 1 && lower == 0 && return [(1, 0, 0)]
  lower == 0 && return [(1, p, 0)]
  upper == q + 1 && return [((-1)^s * binomial(p, s), 0, s) for s in 0:p]
  return [((-1)^s * binomial(p, s), p - s, s) for s in 0:p]
end

function apply_powers(
  monomial::TimeMonomial, upper::Int, lower::Int, a::Int, b::Int, q::Int
)
  result = monomial
  a > 0 && upper <= q && (result = shifted_monomial(result, upper, 0, a))
  b > 0 && lower >= 1 && (result = shifted_monomial(result, lower, 0, b))
  return result
end

function drift(
  kernel::TimeKernel{X},
  powers::Vector{Vector{X}},
  twopi,
  upper::Int,
  lower::Int,
  q::Int,
  R::Int,
) where {X}
  terms = Dict{KernelKey,X}()
  for ((order, monomial), value) in kernel.terms, p in 0:(R - order)
    scale = twopi^p * (1 // factorial(p))
    for (n, Qp) in pairs(powers[p + 1])
      total = order + p + n - 1
      total <= R || break
      iszero(Qp) && continue
      product = Qp * value
      for (c, a, b) in difference_expansion(p, upper, lower, q)
        key = (total, apply_powers(monomial, upper, lower, a, b, q))
        accumulate!(terms, key, (c * scale) * product)
      end
    end
  end
  return TimeKernel{X}(terms)
end

function insert_jump(
  kernel::TimeKernel{X}, jump::Vector{Harmonics{X}}, variable::Int, R::Int
) where {X}
  terms = Dict{KernelKey,X}()
  for ((order, monomial), value) in kernel.terms, (n, harmonics) in pairs(jump)
    total = order + n - 1
    total <= R || continue
    for (k, J) in harmonics
      key = (total, shifted_monomial(monomial, variable, k, 0))
      accumulate!(terms, key, J * value)
    end
  end
  return TimeKernel{X}(terms)
end

function conjugate_kernel(
  kernel::TimeKernel{X}, left::Vector{X}, right::Vector{X}, R::Int
) where {X}
  terms = Dict{KernelKey,X}()
  for ((order, monomial), value) in kernel.terms,
    (n, A) in pairs(left),
    (m, B) in pairs(right)

    total = order + n + m - 2
    total <= R || continue
    accumulate!(terms, (total, monomial), A * value * B)
  end
  return TimeKernel{X}(terms)
end

function simplex_primitive!(
  target::Dict{Tuple{Int,Int},F}, m::Int, p::Int, c::F, inverse_twopi::F
) where {F}
  if iszero(m)
    accumulate!(target, (0, p + 1), c * (1 // (p + 1)))
    return target
  end
  inverse_rate = F(0, 1) * inverse_twopi * (1 // m)
  falling = 1
  for j in 0:p
    accumulate!(target, (m, p - j), c * ((-1)^j * falling) * inverse_rate^(j + 1))
    falling *= p - j
  end
  accumulate!(target, (0, 0), -c * ((-1)^p * factorial(p)) * inverse_rate^(p + 1))
  return target
end

function simplex_integral(
  phases::Vector{Int}, powers::Vector{Int}, inverse_twopi::F
) where {F}
  current = Dict{Tuple{Int,Int},F}((0, 0) => one(F))
  for j in eachindex(phases)
    next = Dict{Tuple{Int,Int},F}()
    for ((m, p), c) in current
      simplex_primitive!(next, m + phases[j], p + powers[j], c, inverse_twopi)
    end
    current = next
  end
  return sum(values(current); init=zero(F))
end

function kernel_gram(monomials::Vector{TimeMonomial}, inverse_twopi::F) where {F}
  n = length(monomials)
  gram = Matrix{F}(undef, n, n)
  for j in 1:n, i in j:n
    (phases_i, powers_i), (phases_j, powers_j) = monomials[i], monomials[j]
    value = simplex_integral(phases_i - phases_j, powers_i + powers_j, inverse_twopi)
    gram[i, j] = value
    gram[j, i] = conj(value)
  end
  return gram
end

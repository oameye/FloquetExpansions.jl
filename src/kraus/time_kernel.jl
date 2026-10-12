const TimeMonomial = Tuple{Vector{Int},Vector{Int}}
const KernelKey = Tuple{Int,TimeMonomial}

struct TimeKernel{X}
  terms::Dict{KernelKey,X}
end

function kernel_identity(identity::X, q::Int) where {X}
  return TimeKernel{X}(Dict{KernelKey,X}((0, (zeros(Int, q), zeros(Int, q))) => identity))
end

function shifted_phases((phases, powers)::TimeMonomial, variable::Int, phase::Int)
  shifted = copy(phases)
  shifted[variable] += phase
  return (shifted, powers)
end

function shifted_powers((phases, powers)::TimeMonomial, variable::Int, power::Int)
  iszero(power) && return (phases, powers)
  shifted = copy(powers)
  shifted[variable] += power
  return (phases, shifted)
end

function segment_expansion(monomial::TimeMonomial, p::Int, upper::Int, lower::Int)
  open_upper = upper <= length(first(monomial))
  lower == 0 && return [(1, open_upper ? shifted_powers(monomial, upper, p) : monomial)]
  return [
    (
      (-1)^s * binomial(p, s),
      shifted_powers(
        open_upper ? shifted_powers(monomial, upper, p - s) : monomial, lower, s
      ),
    ) for s in 0:p
  ]
end

function drift(
  kernel::TimeKernel{X}, powers::Vector{Vector{X}}, twopi, upper::Int, lower::Int, R::Int
) where {X}
  scales = [twopi^p * (1 // factorial(p)) for p in 0:R]
  terms = Dict{KernelKey,X}()
  for ((order, monomial), value) in kernel.terms, p in 0:(R - order)
    expansion = segment_expansion(monomial, p, upper, lower)
    for (n, Qp) in pairs(powers[p + 1])
      total = order + p + n - 1
      total <= R || break
      iszero(Qp) && continue
      product = Qp * value
      for (c, shifted) in expansion
        accumulate!(terms, (total, shifted), (c * scales[p + 1]) * product)
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
      key = (total, shifted_phases(monomial, variable, k))
      accumulate!(terms, key, J * value)
    end
  end
  return TimeKernel{X}(terms)
end

function conjugate_kernel(
  kernel::TimeKernel{X}, left::Vector{X}, right::Vector{X}, R::Int
) where {X}
  terms = Dict{KernelKey,X}()
  for ((order, monomial), value) in kernel.terms, (n, A) in pairs(left)
    order + n - 1 <= R || continue
    Avalue = A * value
    for (m, B) in pairs(right)
      total = order + n + m - 2
      total <= R || continue
      accumulate!(terms, (total, monomial), Avalue * B)
    end
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

struct SimplexIntegrals{F}
  inverse_twopi::F
  values::Dict{TimeMonomial,F}
end

function SimplexIntegrals(inverse_twopi::F) where {F}
  return SimplexIntegrals(inverse_twopi, Dict{TimeMonomial,F}())
end

function gram_entry(integrals::SimplexIntegrals, μ::TimeMonomial, ν::TimeMonomial)
  (phases_μ, powers_μ), (phases_ν, powers_ν) = μ, ν
  return get!(integrals.values, (phases_μ - phases_ν, powers_μ + powers_ν)) do
    return simplex_integral(
      phases_μ - phases_ν, powers_μ + powers_ν, integrals.inverse_twopi
    )
  end
end

const TimeMonomial = Tuple{Vector{Int},Vector{Int}}
const KernelKey = Tuple{Int,Int,TimeMonomial}
const PeriodPolynomial{F} = Dict{Int,F}

struct TimeKernel{X}
  terms::Dict{KernelKey,X}
end

function kernel_identity(identity::X, q::Int) where {X}
  return TimeKernel{X}(
    Dict{KernelKey,X}((0, 0, (zeros(Int, q), zeros(Int, q))) => identity)
  )
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

function raise_upper(monomial::TimeMonomial, upper::Int, power::Int)
  upper <= length(first(monomial)) || return (power, monomial)
  return (0, shifted_powers(monomial, upper, power))
end

function segment_expansion(monomial::TimeMonomial, p::Int, upper::Int, lower::Int)
  lower == 0 && return [(1, raise_upper(monomial, upper, p)...)]
  return map(0:p) do s
    endpoint, raised = raise_upper(monomial, upper, p - s)
    return ((-1)^s * binomial(p, s), endpoint, shifted_powers(raised, lower, s))
  end
end

function drift(
  kernel::TimeKernel{X}, powers::Vector{Vector{X}}, upper::Int, lower::Int, R::Int
) where {X}
  terms = Dict{KernelKey,X}()
  for ((order, period, monomial), value) in kernel.terms, p in 0:(R - order)
    expansion = segment_expansion(monomial, p, upper, lower)
    for (n, Qp) in pairs(powers[p + 1])
      total = order + p + n - 1
      total <= R || break
      iszero(Qp) && continue
      product = Qp * value
      for (c, endpoint, shifted) in expansion
        accumulate!(
          terms, (total, period + endpoint, shifted), (c // factorial(p)) * product
        )
      end
    end
  end
  return TimeKernel{X}(terms)
end

function insert_jump(
  kernel::TimeKernel{X}, jump::Vector{Harmonics{X}}, variable::Int, R::Int
) where {X}
  terms = Dict{KernelKey,X}()
  for ((order, period, monomial), value) in kernel.terms, (n, harmonics) in pairs(jump)
    total = order + n - 1
    total <= R || continue
    for (k, J) in harmonics
      key = (total, period, shifted_phases(monomial, variable, k))
      accumulate!(terms, key, J * value)
    end
  end
  return TimeKernel{X}(terms)
end

function conjugate_kernel(
  kernel::TimeKernel{X}, left::Vector{X}, right::Vector{X}, R::Int
) where {X}
  terms = Dict{KernelKey,X}()
  for ((order, period, monomial), value) in kernel.terms, (n, A) in pairs(left)
    order + n - 1 <= R || continue
    Avalue = A * value
    for (m, B) in pairs(right)
      total = order + n + m - 2
      total <= R || continue
      accumulate!(terms, (total, period, monomial), Avalue * B)
    end
  end
  return TimeKernel{X}(terms)
end

function simplex_primitive!(target::Dict{Tuple{Int,Int},F}, m::Int, p::Int, c::F) where {F}
  if iszero(m)
    accumulate!(target, (0, p + 1), c * (1 // (p + 1)))
    return target
  end
  inverse_rate = F(0, 1) * (1 // m)
  falling = 1
  for j in 0:p
    accumulate!(target, (m, p - j), c * ((-1)^j * falling) * inverse_rate^(j + 1))
    falling *= p - j
  end
  accumulate!(target, (0, 0), -c * ((-1)^p * factorial(p)) * inverse_rate^(p + 1))
  return target
end

function simplex_integral(::Type{F}, (phases, powers)::TimeMonomial) where {F}
  current = Dict{Tuple{Int,Int},F}((0, 0) => one(F))
  for j in eachindex(phases)
    next = Dict{Tuple{Int,Int},F}()
    for ((m, p), c) in current
      simplex_primitive!(next, m + phases[j], p + powers[j], c)
    end
    current = next
  end
  integral = PeriodPolynomial{F}()
  for ((_, p), c) in current
    accumulate!(integral, p, c)
  end
  return integral
end

struct SimplexIntegrals{F}
  values::Dict{TimeMonomial,PeriodPolynomial{F}}
end

SimplexIntegrals{F}() where {F} = SimplexIntegrals(Dict{TimeMonomial,PeriodPolynomial{F}}())

function gram_entry(
  integrals::SimplexIntegrals{F},
  (phases_μ, powers_μ)::TimeMonomial,
  (phases_ν, powers_ν)::TimeMonomial,
) where {F}
  product = (phases_μ - phases_ν, powers_μ + powers_ν)
  return get!(() -> simplex_integral(F, product), integrals.values, product)
end

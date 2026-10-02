struct LoweredSite
  kind::Int
  space_index::Int32
  name_id::Int32
  levels::Int
  ground::Int
  pauli::Bool
end

struct SQALowering{T,R}
  algebra::OperatorAlgebra{T,R}
  sites::Vector{LoweredSite}
  parameters::Dict{Any,Any}
end

function lowered_site_kind(op::SQA.Op)
  (SQA.is_destroy(op) || SQA.is_create(op)) && return BOSON_SITE
  (SQA.is_position(op) || SQA.is_momentum(op)) && return PHASE_SITE
  (SQA.is_transition(op) || SQA.is_pauli(op)) && return LEVEL_SITE
  SQA.is_spin(op) && return SPIN_SITE
  return throw(
    ArgumentError(
      "operator kind $(SQA.optype(op)) is not supported by the native expansion"
    ),
  )
end

function collect_lowered_sites!(sites::Dict{Int32,LoweredSite}, q::SQA.QAdd)
  isempty(q.indices) ||
    throw(ArgumentError("the native expansion does not support symbolic index sums"))
  for (term, _) in q.arguments, op in term.ops
    kind = lowered_site_kind(op)
    levels = if SQA.is_transition(op)
      Int(op.nlev)
    elseif SQA.is_pauli(op)
      2
    else
      0
    end
    ground = SQA.is_transition(op) ? Int(op.g) : 0
    site = LoweredSite(kind, op.space_index, op.name_id, levels, ground, SQA.is_pauli(op))
    previous = get(sites, op.space_index, site)
    (previous.kind, previous.levels, previous.pauli) == (kind, levels, site.pauli) ||
      throw(ArgumentError("one Hilbert subspace carries incompatible operator kinds"))
    sites[op.space_index] = previous
  end
  return sites
end

function SQALowering{T}(
  generator::PeriodicGenerator{Liouvillian}, parameters::AbstractDict
) where {T}
  R = real(T)
  found = Dict{Int32,LoweredSite}()
  for (_, L) in generator.components, (left, right, _) in terms(L)
    collect_lowered_sites!(found, qadd(left))
    collect_lowered_sites!(found, qadd(right))
  end
  sites = sort!(collect(values(found)); by=site -> site.space_index)
  algebra_sites = AlgebraSite{R}[
    if site.kind == LEVEL_SITE
      level_site(fill(one(R), site.levels))
    elseif site.kind == PHASE_SITE
      phase_site(R)
    elseif site.kind == SPIN_SITE
      spin_site(R)
    else
      boson_site(R)
    end for site in sites
  ]
  return SQALowering{T,R}(
    OperatorAlgebra{T}(algebra_sites), sites, Dict{Any,Any}(parameters)
  )
end

function site_position(lowering::SQALowering, op::SQA.Op)
  position = findfirst(site -> site.space_index == op.space_index, lowering.sites)
  position === nothing && throw(ArgumentError("operator acts on an unknown subspace"))
  return position
end

function site_monomial(lowering::SQALowering, s::Int, exponents::Int...)
  monomial = identity_monomial(lowering.algebra)
  monomial[site_slots(lowering.algebra, s)] .= exponents
  return monomial
end

function level_unit(lowering::SQALowering{T}, s::Int, i::Int, j::Int) where {T}
  site = lowering.algebra.sites[s]
  pairs = [
    (site_monomial(lowering, s, m[1], m[2]), c) for
    (m, c) in site_identity_reduction(site, i, j, T)
  ]
  return algebra_operator(lowering.algebra, pairs)
end

# S_x = (S₊ + S₋) / 2, S_y = (S₊ - S₋) / 2i, with S± = S_x ± i S_y
function lower_spin(lowering::SQALowering{T}, s::Int, axis::Int) where {T}
  algebra = lowering.algebra
  axis == 3 && return algebra_operator(algebra, [(site_monomial(lowering, s, 0, 1, 0), 1)])
  raising = site_monomial(lowering, s, 1, 0, 0)
  lowered = site_monomial(lowering, s, 0, 0, 1)
  half = T(1 // 2)
  axis == 1 && return algebra_operator(algebra, [(raising, half), (lowered, half)])
  return algebra_operator(algebra, [(raising, -im * half), (lowered, im * half)])
end

function lower_operator(lowering::SQALowering{T}, op::SQA.Op) where {T}
  s = site_position(lowering, op)
  algebra = lowering.algebra
  SQA.is_destroy(op) &&
    return algebra_operator(algebra, [(site_monomial(lowering, s, 0, 1), 1)])
  SQA.is_create(op) &&
    return algebra_operator(algebra, [(site_monomial(lowering, s, 1, 0), 1)])
  SQA.is_position(op) &&
    return algebra_operator(algebra, [(site_monomial(lowering, s, 1, 0), 1)])
  SQA.is_momentum(op) &&
    return algebra_operator(algebra, [(site_monomial(lowering, s, 0, 1), 1)])
  SQA.is_transition(op) && return level_unit(lowering, s, Int(op.l1), Int(op.l2))
  SQA.is_spin(op) && return lower_spin(lowering, s, Int(op.l1))
  up = level_unit(lowering, s, 1, 2)
  down = level_unit(lowering, s, 2, 1)
  axis = Int(op.l1)
  axis == 1 && return up + down
  axis == 2 && return -im * up + im * down
  return level_unit(lowering, s, 1, 1) - level_unit(lowering, s, 2, 2)
end

exact_real(::Type{R}, x::Integer) where {R} = R(x)
exact_real(::Type{R}, x::Rational) where {R} = R(x)
exact_real(::Type{R}, x::AbstractFloat) where {R} = R(rationalize(BigInt, x))
exact_scalar(::Type{T}, x::Real) where {T} = T(exact_real(real(T), x))
function exact_scalar(::Type{T}, x::Complex) where {T}
  return T(exact_real(real(T), real(x)), exact_real(real(T), imag(x)))
end

function parameter_value(lowering::SQALowering{T}, symbol) where {T}
  key = Symbolics.unwrap(symbol)
  for (k, v) in lowering.parameters
    isequal(Symbolics.unwrap(k), key) && return exact_scalar(T, v)
  end
  return throw(
    ArgumentError("the native expansion needs a numeric value for parameter $symbol")
  )
end

function exact_root(x::Integer, q::Int)
  x < 0 && isodd(q) && return -exact_root(-x, q)
  x < 0 && throw(ArgumentError("an even root of a negative parameter is not real"))
  r = round(BigInt, big(x)^(1 / q))
  for candidate in (r - 1, r, r + 1)
    candidate >= 0 && candidate^q == x && return candidate
  end
  return throw(
    ArgumentError("the root of a parameter is irrational; supply exact jump rates instead")
  )
end

function exact_power(x::T, exponent::Rational) where {T}
  denominator(exponent) == 1 && return x^numerator(exponent)
  iszero(imag(x)) || throw(ArgumentError("fractional powers need real parameters"))
  r = real(x)
  q = denominator(exponent)
  root = exact_root(numerator(r), q) // exact_root(denominator(r), q)
  return T(root)^numerator(exponent)
end

function lower_coefficient(lowering::SQALowering{T}, c::SQA.CNum) where {T}
  tail = c.tail
  tail isa SQA.Native && return exact_scalar(T, c.z)
  tail isa SQA.Poly || throw(
    ArgumentError("the native expansion supports polynomial coefficients only, got $c")
  )
  value = zero(T)
  for monomial in tail.terms
    term = exact_scalar(T, monomial.scalar)
    for (symbol, exponent) in zip(monomial.syms, monomial.exps)
      term *= exact_power(parameter_value(lowering, symbol), exponent)
    end
    value += term
  end
  return value
end

function lower_qadd(lowering::SQALowering{T,R}, q::SQA.QAdd) where {T,R}
  result = AlgebraOperator{T,R}(lowering.algebra, Dict{Monomial,T}())
  for (term, c) in q.arguments
    product = algebra_identity(lowering.algebra)
    for op in term.ops
      product = product * lower_operator(lowering, op)
    end
    result = result + lower_coefficient(lowering, c) * product
  end
  return result
end

function lower_term(lowering::SQALowering{T,R}, term::SQA.QTerm) where {T,R}
  product = algebra_identity(lowering.algebra)
  for op in term.ops
    product = product * lower_operator(lowering, op)
  end
  return product
end

function lowered_monomials(lowering::SQALowering{T,R}, q::SQA.QAdd) where {T,R}
  result = Tuple{AlgebraOperator{T,R},SQA.Coeff}[]
  for (term, c) in q.arguments
    push!(result, (lower_term(lowering, term), c))
  end
  return result
end

function lower_liouvillian(lowering::SQALowering{T,R}, L::Liouvillian) where {T,R}
  result = AlgebraSuperoperator{T,R}(lowering.algebra, Dict{Tuple{Monomial,Monomial},T}())
  for (left, right, c) in terms(L)
    for (X, cl) in lowered_monomials(lowering, qadd(left)),
      (Y, cr) in lowered_monomials(lowering, qadd(right))

      result = result + lower_coefficient(lowering, c * cl * cr) * sandwich(X, Y)
    end
  end
  return result
end

function lower_generator(
  lowering::SQALowering{T,R}, generator::PeriodicGenerator{Liouvillian}
) where {T,R}
  return Dict{Int,AlgebraSuperoperator{T,R}}(
    k => lower_liouvillian(lowering, L) for (k, L) in generator.components
  )
end

function spin_site_operator(site::LoweredSite, a::Int, b::Int, c::Int)
  name = SQA.name_from_id(site.name_id)
  index = Int(site.space_index)
  x = qadd(SQA.Spin(name, 1, index))
  y = qadd(SQA.Spin(name, 2, index))
  z = qadd(SQA.Spin(name, 3, index))
  return (x + im * y)^a * z^b * (x - im * y)^c
end

function site_operator(lowering::SQALowering, s::Int, p::Int, q::Int, r::Int...)
  site = lowering.sites[s]
  name = SQA.name_from_id(site.name_id)
  index = Int(site.space_index)
  if site.kind == BOSON_SITE
    a = SQA.Destroy(name, index)
    return qadd(adjoint(a))^p * qadd(a)^q
  elseif site.kind == PHASE_SITE
    return qadd(SQA.Position(name, index))^p * qadd(SQA.Momentum(name, index))^q
  elseif site.kind == SPIN_SITE
    return spin_site_operator(site, p, q, r[1])
  end
  (p, q) == (0, 0) && return one(SQA.QAdd)
  if site.pauli
    x = qadd(SQA.Pauli(name, 1, index))
    y = qadd(SQA.Pauli(name, 2, index))
    z = qadd(SQA.Pauli(name, 3, index))
    (p, q) == (1, 1) && return (1 + z) / 2
    (p, q) == (1, 2) && return (x + im * y) / 2
    (p, q) == (2, 1) && return (x - im * y) / 2
    return (1 - z) / 2
  end
  return qadd(SQA.Transition(name, p, q, index, SQA.NO_INDEX, site.ground, site.levels))
end

function sqa_scalar(x::Complex{<:Rational})
  try
    return Complex{Rational{Int}}(x)
  catch error
    error isa InexactError || rethrow()
    return ComplexF64(x)
  end
end
sqa_scalar(x::Number) = x

function lift_operator(lowering::SQALowering, X::AlgebraOperator)
  result = zero(SQA.QAdd)
  for (monomial, c) in X.terms
    product = one(SQA.QAdd)
    for s in eachindex(lowering.sites)
      slots = monomial[site_slots(lowering.algebra, s)]
      product = product * site_operator(lowering, s, slots...)
    end
    result = result + sqa_scalar(c) * product
  end
  return SQA.simplify(result)
end

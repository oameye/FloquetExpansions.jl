const MAXIMUM_PARAMETERS = 8
const PARAMETER_NAMES = [Symbol(:p, i) for i in 1:MAXIMUM_PARAMETERS]
const ParameterFraction = Nemo.Generic.FracFieldElem{Nemo.QQMPolyRingElem}

struct ParameterScalar <: Real
  value::ParameterFraction
  ParameterScalar(value::ParameterFraction) = new(value)
end

parameter_ring() = Nemo.QQMPolyRing(Nemo.QQ, PARAMETER_NAMES, :lex, true)

parameter_field() = Nemo.fraction_field(parameter_ring(); cached=true)

function parameter_variable(i::Int)
  return ParameterScalar(parameter_field()(Nemo.gen(parameter_ring(), i)))
end

ParameterScalar(x::Bool) = ParameterScalar(Int(x))
ParameterScalar(x::Integer) = ParameterScalar(parameter_field()(Nemo.QQ(x)))
function ParameterScalar(x::Rational)
  return ParameterScalar(parameter_field()(Nemo.QQ(numerator(x), denominator(x))))
end

Base.promote_rule(::Type{ParameterScalar}, ::Type{<:Integer}) = ParameterScalar
Base.promote_rule(::Type{ParameterScalar}, ::Type{<:Rational}) = ParameterScalar

Base.zero(::Type{ParameterScalar}) = ParameterScalar(0)
Base.one(::Type{ParameterScalar}) = ParameterScalar(1)

function Base.:+(a::ParameterScalar, b::ParameterScalar)
  iszero(a) && return b
  iszero(b) && return a
  return ParameterScalar(a.value + b.value)
end

function Base.:-(a::ParameterScalar, b::ParameterScalar)
  iszero(b) && return a
  iszero(a) && return -b
  return ParameterScalar(a.value - b.value)
end

function Base.:*(a::ParameterScalar, b::ParameterScalar)
  iszero(a) && return a
  iszero(b) && return b
  isone(a) && return b
  isone(b) && return a
  return ParameterScalar(a.value * b.value)
end

Base.:/(a::ParameterScalar, b::ParameterScalar) = ParameterScalar(a.value // b.value)
Base.:-(a::ParameterScalar) = ParameterScalar(-a.value)
Base.inv(a::ParameterScalar) = ParameterScalar(inv(a.value))
Base.:^(a::ParameterScalar, n::Integer) = ParameterScalar(a.value^Int(n))

Base.iszero(a::ParameterScalar) = iszero(a.value)
Base.isone(a::ParameterScalar) = isone(a.value)
Base.:(==)(a::ParameterScalar, b::ParameterScalar) = a.value == b.value
Base.isequal(a::ParameterScalar, b::ParameterScalar) = a.value == b.value
Base.hash(a::ParameterScalar, h::UInt) = hash(a.value, h)
Base.isfinite(::ParameterScalar) = true
Base.isinf(::ParameterScalar) = false
Base.isnan(::ParameterScalar) = false
Base.show(io::IO, a::ParameterScalar) = print(io, a.value)

struct ParameterSignError <: Exception
  value::ParameterScalar
end

function definite_sign(p::Nemo.QQMPolyRingElem, a::ParameterScalar)
  signs = unique(Nemo.sign(Int, c) for c in Nemo.coefficients(p))
  length(signs) == 1 || throw(ParameterSignError(a))
  return only(signs)
end

function parameter_sign(a::ParameterScalar)
  iszero(a) && return 0
  return definite_sign(numerator(a.value), a) * definite_sign(denominator(a.value), a)
end

Base.:<(a::ParameterScalar, b::ParameterScalar) = parameter_sign(b - a) > 0
Base.:<=(a::ParameterScalar, b::ParameterScalar) = parameter_sign(b - a) >= 0
Base.isless(a::ParameterScalar, b::ParameterScalar) = a < b

Base.inv(z::Complex{ParameterScalar}) = conj(z) / abs2(z)
Base.:/(a::Complex{ParameterScalar}, b::Complex{ParameterScalar}) = a * inv(b)
function Base.:^(z::Complex{ParameterScalar}, n::Int)
  n >= 0 && return Base.power_by_squaring(z, n)
  return inv(Base.power_by_squaring(z, -n))
end

struct ParameterCoefficients
  parameters::Vector{Symbolics.Num}
  signs::Vector{Int}
end

function collect_parameters!(found::Vector{Symbolics.Num}, c::SQA.CNum)
  tail = c.tail
  tail isa SQA.Native && return found
  tail isa SQA.Poly || throw(
    ArgumentError("the native expansion supports polynomial coefficients only, got $c")
  )
  for monomial in tail.terms, (symbol, exponent) in zip(monomial.syms, monomial.exps)
    isinteger(exponent) || throw(
      ArgumentError(
        "the native expansion supports integer powers of a parameter only, found " *
        "$symbol^($exponent)",
      ),
    )
    parameter = Symbolics.Num(symbol)
    any(p -> isequal(p, parameter), found) || push!(found, parameter)
  end
  return found
end

function collect_parameters!(found::Vector{Symbolics.Num}, q::SQA.QAdd)
  for (_, c) in q.arguments
    collect_parameters!(found, c)
  end
  return found
end

function native_parameters(generator::PeriodicGenerator{Liouvillian})
  found = Symbolics.Num[]
  for (_, L) in generator.components, (left, right, c) in terms(L)
    collect_parameters!(found, c)
    collect_parameters!(found, qadd(left))
    collect_parameters!(found, qadd(right))
  end
  return sort!(found; by=string)
end

function native_coefficients(
  generator::PeriodicGenerator{Liouvillian}, negative::Vector{Symbolics.Num}
)
  parameters = native_parameters(generator)
  for p in negative
    isequal(p, generator.wd) && throw(
      ArgumentError("the drive frequency is positive and cannot be declared negative")
    )
    any(q -> isequal(q, p), parameters) || throw(
      ArgumentError("the parameter $p declared negative does not occur in the generator")
    )
  end
  isempty(parameters) && return ExactCoefficients()
  length(parameters) <= MAXIMUM_PARAMETERS || throw(
    ArgumentError(
      "the GKSL normal form supports at most $MAXIMUM_PARAMETERS symbolic parameters, " *
      "found $(length(parameters)); substitute exact values for some of them",
    ),
  )
  signs = [any(q -> isequal(q, p), negative) ? -1 : 1 for p in parameters]
  return ParameterCoefficients(parameters, signs)
end

function native_expansion_data(
  ::ExactCoefficients,
  algorithm::Union{BlochFeshbach,HoriDeprit},
  generator::PeriodicGenerator{Liouvillian},
  N::Int,
  gauge_degree::Int,
)
  return native_expansion_data(algorithm, generator, N, gauge_degree)
end

function native_expansion_data(
  coefficients::ParameterCoefficients,
  algorithm::Union{BlochFeshbach,HoriDeprit},
  generator::PeriodicGenerator{Liouvillian},
  N::Int,
  gauge_degree::Int,
)
  T = Complex{ParameterScalar}
  values = T[s * parameter_variable(i) for (i, s) in pairs(coefficients.signs)]
  lowering = SQALowering{T}(generator, coefficients.parameters, values)
  try
    return native_expansion_data(lowering, algorithm, generator, N, gauge_degree)
  catch error
    error isa ParameterSignError || rethrow()
    value = parameter_coefficient(coefficients, error.value)
    throw(
      ArgumentError(
        "the sign of $value is not fixed by the declared parameter signs; pass `negative` " *
        "to GKSLNormalForm or substitute exact values for some parameters",
      ),
    )
  end
end

function parameter_term(coefficients::ParameterCoefficients, i::Int)
  return convert(SQA.CNum, coefficients.signs[i] * coefficients.parameters[i])
end

function parameter_polynomial(coefficients::ParameterCoefficients, p::Nemo.QQMPolyRingElem)
  result = convert(SQA.CNum, 0)
  for (q, exponents) in zip(Nemo.coefficients(p), Nemo.exponent_vectors(p))
    term = convert(SQA.CNum, Rational{BigInt}(q))
    for (i, e) in pairs(exponents)
      iszero(e) || (term *= parameter_term(coefficients, i)^e)
    end
    result += term
  end
  return result
end

function parameter_coefficient(coefficients::ParameterCoefficients, x::ParameterScalar)
  top = parameter_polynomial(coefficients, numerator(x.value))
  bottom = denominator(x.value)
  Nemo.is_constant(bottom) || return top / parameter_polynomial(coefficients, bottom)
  scale = inv(Rational{BigInt}(Nemo.constant_coefficient(bottom)))
  return top * convert(SQA.CNum, scale)
end

function lift_coefficient(coefficients::ParameterCoefficients, c::Complex{ParameterScalar})
  real_part = parameter_coefficient(coefficients, real(c))
  iszero(imag(c)) && return real_part
  imaginary_part = parameter_coefficient(coefficients, imag(c))
  return real_part + convert(SQA.CNum, im) * imaginary_part
end

native_positivity_conditions(::ExactCoefficients, jumps, wd::Symbolics.Num) = SQA.CNum[]

function native_positivity_conditions(
  coefficients::ParameterCoefficients, jumps, wd::Symbolics.Num
)
  conditions = SQA.CNum[
    parameter_term(coefficients, i) for i in eachindex(coefficients.parameters)
  ]
  any(p -> isequal(p, wd), coefficients.parameters) ||
    push!(conditions, convert(SQA.CNum, wd))
  for jump in jumps
    any(c -> isequal(c, jump.rate), conditions) || push!(conditions, jump.rate)
  end
  return conditions
end

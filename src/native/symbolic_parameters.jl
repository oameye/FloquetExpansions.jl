function collect_native_parameters!(found::Vector{Any}, c::SQA.CNum)
  tail = c.tail
  tail isa SQA.Native && return found
  tail isa SQA.Poly || throw(
    ArgumentError("the native expansion supports polynomial coefficients only, got $c")
  )
  for monomial in tail.terms, symbol in monomial.syms
    any(s -> isequal(s, symbol)::Bool, found) || push!(found, symbol)
  end
  return found
end

function native_parameters(generator::PeriodicGenerator{Liouvillian})
  found = Any[]
  for (_, L) in generator.components, (_, _, c) in terms(L)
    collect_native_parameters!(found, c)
  end
  return found
end

# Generic structure points: ratios of distinct primes, away from the small-integer
# relations at which resonances and accidental cancellations occur.
const STRUCTURE_PRIMES = (101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163)

function structure_values(m::Int, k::Int)
  primes = STRUCTURE_PRIMES
  n = length(primes)
  return Rational{BigInt}[
    primes[mod1(2i + k, n)] // primes[mod1(3i + 2k + 1, n)] for i in 1:m
  ]
end

function parameter_point(parameters::Vector{Any}, values::Vector{Rational{BigInt}})
  return Dict{Any,Any}(p => x for (p, x) in zip(parameters, values))
end

function structure_point(parameters::Vector{Any}, k::Int)
  return parameter_point(parameters, structure_values(length(parameters), k))
end

const NativeOutputKey = Tuple{Symbol,Int,Int,Monomial,Monomial}

function record_superoperator!(outputs::Dict, kind::Symbol, n::Int, k::Int, S)
  for ((X, Y), c) in S.terms
    iszero(c) || (outputs[(kind, n, k, X, Y)] = c)
  end
  return outputs
end

function record_operator!(outputs::Dict, kind::Symbol, n::Int, k::Int, X::AlgebraOperator)
  for (M, c) in X.terms
    iszero(c) || (outputs[(kind, n, k, M, Monomial())] = c)
  end
  return outputs
end

function record_matrix!(outputs::Dict, kind::Symbol, n::Int, C::AbstractMatrix)
  for index in CartesianIndices(C)
    μ, ν = Tuple(index)
    iszero(C[index]) || (outputs[(kind, n, 0, [μ], [ν])] = C[index])
  end
  return outputs
end

function record_recurrence!(outputs::Dict, data::NativeExpansionData)
  for (n, E) in pairs(data.recurrence.E)
    record_superoperator!(outputs, :effective, n, 0, E)
    record_matrix!(outputs, :retained, n, native_kossakowski(data.representation, E))
  end
  for (n, H) in pairs(data.recurrence.H)
    record_operator!(outputs, :hamiltonian, n, 0, H)
  end
  return outputs
end

function record_micromotion!(outputs::Dict, data::NativeExpansionData, algorithm, N::Int)
  for (n, harmonics) in pairs(native_micromotion_harmonics(data, algorithm, N)),
    (k, X) in harmonics

    record_superoperator!(outputs, :micromotion, n, k, X)
  end
  return outputs
end

function record_channels!(outputs::Dict, data::NativeExpansionData, N::Int)
  for (j, channel) in pairs(data.recurrence.channels)
    channel.onset <= N || continue
    outputs[(:weight, j, 0, Monomial(), Monomial())] = channel.weight
    for (order, coefficients) in retained_amplitudes(channel, N),
      (μ, c) in pairs(coefficients)

      iszero(c) || (outputs[(:amplitude, j, order, [μ], Monomial())] = c)
    end
  end
  return outputs
end

function native_outputs(data::NativeExpansionData{T}, algorithm, N::Int) where {T}
  outputs = Dict{NativeOutputKey,T}()
  record_recurrence!(outputs, data)
  record_micromotion!(outputs, data, algorithm, N)
  record_channels!(outputs, data, N)
  graded = native_finite_kossakowski(data, N)
  for (order, C) in graded
    record_matrix!(outputs, :graded, order, C)
  end
  for (order, G) in native_finite_superoperators(data, graded, N)
    record_superoperator!(outputs, :finite, order, 0, G)
  end
  return outputs
end

function static_gauge_orders(data::NativeExpansionData)
  recurrence = data.recurrence
  return [
    n for n in eachindex(recurrence.S) if
    !iszero(recurrence.S[n]) || haskey(recurrence.virtual, n)
  ]
end

const NativeStructure = Tuple{
  Vector{Monomial},Vector{Int},Vector{Int},Vector{Int},Set{NativeOutputKey}
}

function native_structure(data::NativeExpansionData, outputs::Dict)::NativeStructure
  return (
    copy(data.representation.frame),
    [channel.onset for channel in data.recurrence.channels],
    sort!(collect(keys(data.recurrence.virtual))),
    static_gauge_orders(data),
    Set{NativeOutputKey}(keys(outputs)),
  )
end

struct NativeSamples{A<:ExpansionAlgorithm}
  generator::PeriodicGenerator{Liouvillian}
  algorithm::A
  order::Int
  parameters::Vector{Any}
  structure::NativeStructure
end

function structure_mismatch_error(values::Vector{Rational{BigInt}})
  return ArgumentError(
    "the symbolic GKSL normal form changes structure (frame, channels, static gauges, or " *
    "nonzero coefficients) at the parameter point $(join(string.(values), ", ")) of the " *
    "sign region of the reference point. The result is piecewise in the parameters; " *
    "substitute exact values instead.",
  )
end

function (samples::NativeSamples)(values::Vector{Rational{BigInt}})
  point = parameter_point(samples.parameters, values)
  data = native_expansion_data(
    samples.algorithm, samples.generator, samples.order, point, 3
  )
  outputs = native_outputs(data, samples.algorithm, samples.order)
  isequal(native_structure(data, outputs), samples.structure) ||
    throw(structure_mismatch_error(values))
  return Dict{NativeOutputKey,ExactField}(k => ExactField(v) for (k, v) in outputs)
end

function parameter_coefficient(parameters::Vector{Any}, exponents::Vector{Int})
  value = convert(SQA.CNum, 1)
  for (p, e) in zip(parameters, exponents)
    e == 0 || (value = value * convert(SQA.CNum, Symbolics.Num(p))^e)
  end
  return value
end

function sqa_polynomial(
  parameters::Vector{Any}, terms::Vector{Tuple{Vector{Int},ExactField}}
)
  value = convert(SQA.CNum, 0)
  for (exponents, c) in terms
    value =
      value +
      convert(SQA.CNum, sqa_scalar(c)) * parameter_coefficient(parameters, exponents)
  end
  return value
end

function sqa_rational(parameters::Vector{Any}, f::RationalFunction)
  numerator = sqa_polynomial(parameters, f.numerator)
  length(f.denominator) == 1 &&
    iszero(sum(f.denominator[1][1])) &&
    return numerator * convert(SQA.CNum, sqa_scalar(inv(f.denominator[1][2])))
  return numerator / sqa_polynomial(parameters, f.denominator)
end

struct ReconstructedCoefficients{D<:NativeExpansionData,T}
  reference::D
  point::Dict{Any,Any}
  values::Dict{T,SQA.CNum}
  conditions::Vector{SQA.CNum}
end

function (coefficients::ReconstructedCoefficients)(c)
  iszero(c) && return convert(SQA.CNum, 0)
  return get(coefficients.values, c) do
    return throw(
      ArgumentError("a coefficient of the symbolic GKSL normal form was not reconstructed")
    )
  end
end

function native_positivity_conditions(coefficients::ReconstructedCoefficients)
  return copy(coefficients.conditions)
end

function native_factorization(
  ::NativeExpansionData, coefficients::ReconstructedCoefficients
)
  return coefficients
end

function coefficient_values(
  reference::Dict{NativeOutputKey,T}, fitted::Dict, parameters::Vector{Any}
) where {T}
  functions = Dict{T,RationalFunction}()
  values = Dict{T,SQA.CNum}()
  for (key, c) in reference
    f = fitted[key]
    get(functions, c, f) == f || throw(
      ArgumentError(
        "two distinct coefficients of the symbolic GKSL normal form coincide at the " *
        "reference point; substitute exact values instead.",
      ),
    )
    haskey(values, c) && continue
    functions[c] = f
    values[c] = sqa_rational(parameters, f)
  end
  return values
end

function region_conditions(parameters::Vector{Any}, signs::Vector{Int})
  return SQA.CNum[
    convert(SQA.CNum, s) * convert(SQA.CNum, Symbolics.Num(p)) for
    (p, s) in zip(parameters, signs)
  ]
end

function rate_conditions(data::NativeExpansionData, values::Dict, N::Int)
  return SQA.CNum[
    values[channel.weight] for channel in data.recurrence.channels if channel.onset <= N
  ]
end

function reconstructed_coefficients(
  reference::NativeExpansionData,
  generator::PeriodicGenerator{Liouvillian},
  algorithm::ExpansionAlgorithm,
  parameters::Vector{Any},
)
  N = length(reference.recurrence.E) - 1
  point = structure_values(length(parameters), 1)
  signs = Int.(sign.(point))
  outputs = native_outputs(reference, algorithm, N)
  structure = native_structure(reference, outputs)
  samples = NativeSamples(generator, algorithm, N, parameters, structure)
  sampler = ExactSampler(samples, collect(structure[5]), signs)
  fitted = reconstruct_rational_functions(sampler)
  verify_reconstruction(
    fitted,
    point,
    Dict{NativeOutputKey,ExactField}(k => ExactField(v) for (k, v) in outputs),
  )
  values = coefficient_values(outputs, fitted, parameters)
  conditions = vcat(
    region_conditions(parameters, signs), rate_conditions(reference, values, N)
  )
  return ReconstructedCoefficients(
    reference, parameter_point(parameters, point), values, conditions
  )
end

function symbolic_gksl_normal_form(
  generator::PeriodicGenerator{Liouvillian},
  gauge::GKSLNormalForm,
  order::Int,
  provenance::R,
  parameters::Vector{Any},
) where {R<:FloquetProvenance}
  N = order - 1
  reference = native_expansion_data(
    gauge.algorithm, generator, N, structure_point(parameters, 1), 3
  )
  if isempty(static_gauge_orders(reference))
    second = native_expansion_data(
      gauge.algorithm, generator, N, structure_point(parameters, 2), 3
    )
    isempty(static_gauge_orders(second)) &&
      return gram_normal_form(generator, gauge, order, provenance)
  end
  coefficients = reconstructed_coefficients(
    reference, generator, gauge.algorithm, parameters
  )
  return native_floquet_expansion(reference, generator, gauge, provenance, coefficients)
end

function gram_normal_form(
  generator::PeriodicGenerator{Liouvillian},
  gauge::GKSLNormalForm,
  order::Int,
  provenance::R,
) where {R<:FloquetProvenance}
  reference = floquet_expansion_impl(
    generator, VanVleck(gauge.algorithm), order, provenance
  )
  completed = positive_completion(reference, Gram())
  return FloquetExpansion(
    generator,
    getfield(reference, :micromotion_components),
    getfield(reference, :effective_components),
    gauge,
    order,
    getfield(completed, :completion),
    provenance,
  )
end

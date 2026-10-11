const ExactField = Complex{ExactReal}
const MAXIMUM_PARAMETERS = 6

function collect_parameters!(found::Vector{Symbolics.Num}, c::SQA.CNum)
  tail = c.tail
  tail isa SQA.Native && return found
  tail isa SQA.Poly || throw(
    ArgumentError("the native expansion supports polynomial coefficients only, got $c")
  )
  for monomial in tail.terms, (symbol, exponent) in zip(monomial.syms, monomial.exps)
    isinteger(exponent) || throw(
      ArgumentError(
        "the native expansion supports integer powers of a parameter only, " *
        "found $symbol^($exponent)",
      ),
    )
    parameter = Symbolics.Num(symbol)
    any(p -> isequal(p, parameter), found) || push!(found, parameter)
  end
  return found
end

function native_parameters(generator::PeriodicGenerator{Liouvillian})
  found = Symbolics.Num[]
  for (_, L) in generator.components, (_, _, c) in terms(L)
    collect_parameters!(found, c)
  end
  return sort!(found; by=string)
end

function parameter_signs(
  parameters::Vector{Symbolics.Num}, negative::Vector{Symbolics.Num}, wd::Symbolics.Num
)
  length(parameters) <= MAXIMUM_PARAMETERS || throw(
    ArgumentError(
      "the symbolic GKSL normal form supports at most $MAXIMUM_PARAMETERS parameters, " *
      "found $(length(parameters)); substitute exact values for some of them",
    ),
  )
  for p in negative
    isequal(p, wd) && throw(ArgumentError("the drive frequency is positive"))
    any(q -> isequal(q, p), parameters) || throw(
      ArgumentError("the parameter $p declared negative does not occur in the generator")
    )
  end
  return [any(q -> isequal(q, p), negative) ? -1 : 1 for p in parameters]
end

function point_scalar(c::SQA.CNum, parameters::Vector{Symbolics.Num}, point::ExactPoint)
  tail = c.tail
  tail isa SQA.Native && return exact_scalar(ExactField, SQA.native_scalar(c))
  total = zero(ExactField)
  for monomial in tail.terms
    term = exact_scalar(ExactField, monomial.scalar)
    for (symbol, exponent) in zip(monomial.syms, monomial.exps)
      i = findfirst(p -> isequal(Symbolics.unwrap(p), symbol), parameters)
      term *= point[i]^Int(exponent)
    end
    total += term
  end
  return total
end

function numeric_generator(
  generator::PeriodicGenerator{Liouvillian},
  parameters::Vector{Symbolics.Num},
  point::ExactPoint,
)
  components = Dict{Int,Liouvillian}()
  for (k, L) in generator.components
    numeric = Liouvillian(LiouvillianTerms())
    for (left, right, c) in terms(L)
      add_term!(numeric, left, right, convert(SQA.CNum, point_scalar(c, parameters, point)))
    end
    components[k] = numeric
  end
  return PeriodicGenerator(components, generator.wd, zero(Liouvillian))
end

function record!(outputs::Dict, label::Tuple, entries)
  for (key, c) in entries
    iszero(c) || (outputs[(label..., key)] = c)
  end
  return outputs
end

function matrix_entries(C::AbstractMatrix)
  return (Tuple(index) => C[index] for index in CartesianIndices(C))
end

function native_outputs(data::NativeExpansionData{T}, algorithm, N::Int) where {T}
  outputs = Dict{Tuple,T}()
  representation = data.representation
  for (n, E) in pairs(data.recurrence.E)
    record!(outputs, (:effective, n), E.terms)
    record!(outputs, (:retained, n), matrix_entries(native_kossakowski(representation, E)))
  end
  for (n, H) in pairs(data.recurrence.H)
    record!(outputs, (:hamiltonian, n), H.terms)
  end
  for (n, harmonics) in pairs(native_micromotion_harmonics(data, algorithm, N)),
    (k, X) in harmonics

    record!(outputs, (:micromotion, n, k), X.terms)
  end
  for (j, channel) in pairs(data.recurrence.channels)
    channel.onset <= N || continue
    record!(outputs, (:weight, j), [(0, channel.weight)])
    for (k, amplitude) in retained_amplitudes(channel, N)
      record!(outputs, (:amplitude, j, k), pairs(amplitude))
    end
  end
  graded = native_finite_kossakowski(data, N)
  for (order, C) in graded
    record!(outputs, (:graded, order), matrix_entries(C))
  end
  for (order, G) in native_finite_superoperators(data, graded, N)
    record!(outputs, (:finite, order), G.terms)
  end
  return outputs
end

function native_structure(data::NativeExpansionData)
  recurrence = data.recurrence
  return (
    copy(data.representation.frame),
    [channel.onset for channel in recurrence.channels],
    [n for n in eachindex(recurrence.S) if !iszero(recurrence.S[n])],
  )
end

function real_outputs(outputs::Dict)
  parts = Dict{Tuple,ExactReal}()
  for (key, c) in outputs
    iszero(real(c)) || (parts[(key, 0)] = ExactReal(real(c)))
    iszero(imag(c)) || (parts[(key, 1)] = ExactReal(imag(c)))
  end
  return parts
end

struct PointSampler{A<:ExpansionAlgorithm}
  algorithm::A
  generator::PeriodicGenerator{Liouvillian}
  parameters::Vector{Symbolics.Num}
  order::Int
  structure::Tuple{Vector{Monomial},Vector{Int},Vector{Int}}
  keys::Set{Tuple}
end

function structure_error(point::ExactPoint, reason::String)
  return ArgumentError(
    "the symbolic GKSL normal form changes structure ($reason) at the parameter point " *
    "$(join(string.(point), ", ")) of the declared sign region. The result is piecewise " *
    "in the parameters; substitute exact values for some of them instead.",
  )
end

function point_data(algorithm, numeric::PeriodicGenerator{Liouvillian}, N::Int)
  return try
    data = native_expansion_data(Complex{Rational{Int128}}, algorithm, numeric, N, 3)
    data, native_outputs(data, algorithm, N)
  catch error
    error isa Union{OverflowError,InexactError} || rethrow()
    data = native_expansion_data(Complex{Rational{BigInt}}, algorithm, numeric, N, 3)
    data, native_outputs(data, algorithm, N)
  end
end

function (sampler::PointSampler)(point::ExactPoint)
  numeric = numeric_generator(sampler.generator, sampler.parameters, point)
  data, raw = try
    point_data(sampler.algorithm, numeric, sampler.order)
  catch error
    error isa NativePositivityError || rethrow()
    throw(structure_error(point, "no polynomial static gauge exists"))
  end
  native_structure(data) == sampler.structure ||
    throw(structure_error(point, "frame, channel onsets or static gauge orders differ"))
  outputs = real_outputs(raw)
  issubset(keys(outputs), sampler.keys) ||
    throw(structure_error(point, "a coefficient is nonzero that vanishes at the reference"))
  return outputs
end

function sqa_polynomial(parameters::Vector{Symbolics.Num}, p::SparsePolynomial)
  value = convert(SQA.CNum, 0)
  for (exponents, c) in p
    term = convert(SQA.CNum, c)
    for (parameter, e) in zip(parameters, exponents)
      iszero(e) || (term *= convert(SQA.CNum, parameter)^e)
    end
    value += term
  end
  return value
end

struct ReconstructedCoefficients
  values::Dict{ExactField,SQA.CNum}
  conditions::Vector{SQA.CNum}
  certificate::NamedTuple{
    (:parameters, :signs, :points),
    Tuple{Vector{Symbolics.Num},Vector{Int},Vector{ExactPoint}},
  }
end

function (coefficients::ReconstructedCoefficients)(c)
  iszero(c) && return convert(SQA.CNum, 0)
  return coefficients.values[ExactField(c)]
end

native_positivity_conditions(c::ReconstructedCoefficients) = copy(c.conditions)
native_factorization(::NativeExpansionData, c::ReconstructedCoefficients) = c.certificate

function coefficient_values(
  parameters::Vector{Symbolics.Num}, reference::Dict{Tuple,<:Complex}, functions::Dict
)
  found = Dict{ExactField,Tuple{Any,Any}}()
  values = Dict{ExactField,SQA.CNum}()
  for (key, c) in reference
    pair = (get(functions, (key, 0), nothing), get(functions, (key, 1), nothing))
    get!(found, ExactField(c), pair) == pair || throw(
      ArgumentError(
        "two distinct coefficients of the symbolic GKSL normal form coincide at the " *
        "reference point; substitute exact values for some parameters instead.",
      ),
    )
    real_part, imaginary_part = (
      f === nothing ? convert(SQA.CNum, 0) : sqa_polynomial(parameters, f) for f in pair
    )
    values[ExactField(c)] = real_part + convert(SQA.CNum, im) * imaginary_part
  end
  return values
end

function region_conditions(
  parameters::Vector{Symbolics.Num}, signs::Vector{Int}, wd::Symbolics.Num
)
  conditions = SQA.CNum[
    convert(SQA.CNum, s) * convert(SQA.CNum, p) for (p, s) in zip(parameters, signs)
  ]
  any(p -> isequal(p, wd), parameters) || push!(conditions, convert(SQA.CNum, wd))
  return conditions
end

function channel_rates(data::NativeExpansionData, values::Dict, wd::Symbolics.Num, N::Int)
  return SQA.CNum[
    values[ExactField(channel.weight)] *
    convert(SQA.CNum, inverse_drive_power(wd, channel.onset)) for
    channel in data.recurrence.channels if channel.onset <= N
  ]
end

function reference_data(algorithm, numeric::PeriodicGenerator{Liouvillian}, order::Int)
  return try
    point_data(algorithm, numeric, order - 1)
  catch error
    error isa NativePositivityError || rethrow()
    throw(
      ArgumentError(
        "order $order of the GKSL normal form has no polynomial static gauge at the " *
        "reference point of the declared sign region; lower the requested order or " *
        "change the region",
      ),
    )
  end
end

function symbolic_gksl_normal_form(
  generator::PeriodicGenerator{Liouvillian},
  gauge::GKSLNormalForm,
  order::Int,
  provenance::R,
  parameters::Vector{Symbolics.Num},
) where {R<:FloquetProvenance}
  wd = generator.wd
  signs = parameter_signs(parameters, gauge.negative, wd)
  point = signs .* reference_magnitudes(length(parameters))
  N = order - 1
  data, outputs = reference_data(
    gauge.algorithm, numeric_generator(generator, parameters, point), order
  )
  parts = real_outputs(outputs)
  sampler = PointSampler(
    gauge.algorithm,
    generator,
    parameters,
    N,
    native_structure(data),
    Set{Tuple}(keys(parts)),
  )
  functions, verified = reconstruct_polynomials(
    sampler, collect(keys(parts)), signs, gauge.verification
  )
  values = coefficient_values(parameters, outputs, functions)
  coefficients = ReconstructedCoefficients(
    values,
    vcat(region_conditions(parameters, signs, wd), channel_rates(data, values, wd, N)),
    (; parameters, signs, points=verified),
  )
  return native_floquet_expansion(data, generator, gauge, provenance, coefficients)
end

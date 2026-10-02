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

function structure_point(parameters::Vector{Any}, k::Int)
  primes = STRUCTURE_PRIMES
  n = length(primes)
  return Dict{Any,Any}(
    p => primes[mod1(2i + k, n)] // primes[mod1(3i + 2k + 1, n)] for
    (i, p) in pairs(parameters)
  )
end

function symbolic_static_gauge_error(parameters::Vector{Any}, order::Int)
  return ArgumentError(
    "the GKSL normal form of this model needs a static gauge at order $order. The exact " *
    "static gauge is computed for numeric parameters only; substitute exact values for " *
    "$(join(string.(parameters), ", ")).",
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
  for k in 1:2
    point = structure_point(parameters, k)
    data = native_expansion_data(gauge.algorithm, generator, N, point, 3)
    n = findfirst(
      n -> !iszero(data.recurrence.S[n]) || haskey(data.recurrence.virtual, n),
      eachindex(data.recurrence.S),
    )
    n === nothing || throw(symbolic_static_gauge_error(parameters, n))
  end
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

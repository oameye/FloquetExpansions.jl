struct CPConfiguredAlgorithm{A<:ExpansionAlgorithm,CP} <: ExpansionAlgorithm
  algorithm::A
end

function cp_configured_algorithm(algorithm::A, ::Val{false}) where {A<:ExpansionAlgorithm}
  return algorithm
end
function cp_configured_algorithm(algorithm::A, ::Val{true}) where {A<:ExpansionAlgorithm}
  return CPConfiguredAlgorithm{A,true}(algorithm)
end
function cp_configured_algorithm(::A, ::Val{CP}) where {A<:ExpansionAlgorithm,CP}
  return throw(ArgumentError("complete_positive must be Val(true) or Val(false)"))
end

"""
    HoriDeprit(; complete_positive=Val(true))

Select the Hori–Deprit Lie-transform expansion algorithm. For Liouvillian input,
`complete_positive=Val(true)` applies the package's graded [`Gram`](@ref) positive completion to
the canonical Hori–Deprit Van Vleck expansion. `Val(false)` returns the traditional raw canonical
series exactly. Hamiltonian expansions are unchanged by the policy.
"""
function HoriDeprit(; complete_positive::Val{CP}=Val(true)) where {CP}
  return cp_configured_algorithm(HoriDeprit(Val(:raw)), complete_positive)
end

"""
    BlochFeshbach(; complete_positive=Val(true))

Select the Bloch/Feshbach projection-recurrence expansion algorithm. For Liouvillian input,
`complete_positive=Val(true)` applies the package's graded [`Gram`](@ref) positive completion to
the canonical Bloch/Feshbach Van Vleck expansion. `Val(false)` returns the traditional raw
canonical series exactly. Hamiltonian expansions are unchanged by the policy.
"""
function BlochFeshbach(; complete_positive::Val{CP}=Val(true)) where {CP}
  return cp_configured_algorithm(BlochFeshbach(Val(:raw)), complete_positive)
end

function Base.show(io::IO, ::CPConfiguredAlgorithm{A,true}) where {A<:ExpansionAlgorithm}
  return print(io, nameof(A), "()")
end
function Base.show(
  io::IO, ::MIME"text/plain", algorithm::CPConfiguredAlgorithm{A,true}
) where {A}
  return show(io, algorithm)
end

function Base.show(io::IO, gauge::VanVleck{<:CPConfiguredAlgorithm})
  print(io, "VanVleck(algorithm=")
  show(io, gauge.algorithm)
  return print(io, ")")
end
function Base.show(io::IO, ::MIME"text/plain", gauge::VanVleck{<:CPConfiguredAlgorithm})
  return show(io, gauge)
end

function rewrap_algorithm_expansion(
  expansion::FloquetExpansion, gauge::G, provenance::R
) where {G<:Gauge,R<:FloquetProvenance}
  return FloquetExpansion(
    expansion.generator,
    expansion.kick_components,
    expansion.effective_components,
    gauge,
    expansion.order,
    getfield(expansion, :completion),
    provenance,
  )
end

function raw_algorithm_gauge(gauge::VanVleck{<:CPConfiguredAlgorithm})
  return VanVleck(gauge.algorithm.algorithm)
end

function floquet_expansion_impl(
  generator::PeriodicGenerator{SQA.QAdd},
  gauge::VanVleck{CPConfiguredAlgorithm{A,true}},
  order::Int,
  provenance::R,
) where {A<:ExpansionAlgorithm,R<:FloquetProvenance}
  expansion = floquet_expansion_impl(
    generator, raw_algorithm_gauge(gauge), order, provenance
  )
  return rewrap_algorithm_expansion(expansion, gauge, provenance)
end

function floquet_expansion_impl(
  generator::PeriodicGenerator{Liouvillian},
  gauge::VanVleck{CPConfiguredAlgorithm{A,true}},
  order::Int,
  provenance::R,
) where {A<:ExpansionAlgorithm,R<:FloquetProvenance}
  raw = floquet_expansion_impl(generator, raw_algorithm_gauge(gauge), order, provenance)
  configured = rewrap_algorithm_expansion(raw, gauge, provenance)
  return positive_completion(configured, Gram())
end

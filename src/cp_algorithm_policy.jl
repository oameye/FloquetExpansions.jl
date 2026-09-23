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

function HoriDeprit(; complete_positive::Val{CP}=Val(false)) where {CP}
  return cp_configured_algorithm(HoriDeprit(Val(:raw)), complete_positive)
end

function BlochFeshbach(; complete_positive::Val{CP}=Val(false)) where {CP}
  return cp_configured_algorithm(BlochFeshbach(Val(:raw)), complete_positive)
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

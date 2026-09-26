bloch_conventions(::PeriodicGenerator{SQA.QAdd}) = BlochConventions(*, 1, im)
bloch_conventions(::PeriodicGenerator{Liouvillian}) = BlochConventions(compose, im, 1)

function floquet_expansion_impl(
  generator::PeriodicGenerator{SQA.QAdd},
  gauge::VanVleck{BlochFeshbach},
  order::Int,
  provenance::R,
) where {R<:FloquetProvenance}
  return bloch_feshbach_expansion(
    generator, gauge, order, provenance, bloch_conventions(generator)
  )
end

function floquet_expansion_impl(
  generator::PeriodicGenerator{Liouvillian},
  gauge::VanVleck{BlochFeshbach},
  order::Int,
  provenance::R,
) where {R<:FloquetProvenance}
  return bloch_feshbach_expansion(
    generator, gauge, order, provenance, bloch_conventions(generator)
  )
end

function bloch_feshbach_zero_expansion(
  generator::P, gauge::G, order::Int, provenance::R
) where {T,P<:PeriodicGenerator{T},G<:VanVleck{BlochFeshbach},R<:FloquetProvenance}
  kick_components = P[zero(generator) for _ in 1:(order - 1)]
  effective_components = T[generator.zero_component for _ in 1:order]
  return FloquetExpansion(
    generator,
    kick_components,
    effective_components,
    gauge,
    order,
    Uncompleted(),
    provenance,
  )
end

function bloch_feshbach_leading_expansion(
  generator::P, gauge::G, provenance::R
) where {T,P<:PeriodicGenerator{T},G<:VanVleck{BlochFeshbach},R<:FloquetProvenance}
  effective = SQA.simplify(time_average(generator))::T
  return FloquetExpansion(generator, P[], T[effective], gauge, 1, Uncompleted(), provenance)
end

function bloch_feshbach_expansion(
  generator::P, gauge::G, order::Int, provenance::R, conventions::BlochConventions
) where {T,P<:PeriodicGenerator{T},G<:VanVleck{BlochFeshbach},R<:FloquetProvenance}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  generator isa PeriodicGenerator{SQA.QAdd} && require_hermitian_drive(generator)
  order == 1 && return bloch_feshbach_leading_expansion(generator, gauge, provenance)
  isempty(generator) &&
    return bloch_feshbach_zero_expansion(generator, gauge, order, provenance)

  (; components, zero_component, wd) = generator
  wave_plan = compile_wave_operator_plan(collect(keys(generator)), order, 0)
  wave_operator = evaluate_wave_operator(wave_plan, components, zero_component, conventions)
  normalization = normalize_to_van_vleck(
    compile_van_vleck_normalization_plan(wave_plan),
    wave_operator,
    components,
    zero_component,
    conventions,
  )

  kick_components = P[
    SQA.simplify(
      conventions.kick_phase * periodic_generator(series, wd, zero_component)
    )::P for series in normalization.connected_log
  ]
  effective_components = T[
    SQA.simplify(component)::T for component in normalization.effective
  ]

  return FloquetExpansion(
    generator,
    kick_components,
    effective_components,
    gauge,
    order,
    Uncompleted(),
    provenance,
  )
end

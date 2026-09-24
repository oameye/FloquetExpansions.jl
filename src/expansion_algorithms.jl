bloch_hamiltonian_product(left::SQA.QAdd, right::SQA.QAdd) = left * right
bloch_liouvillian_product(left::Liouvillian, right::Liouvillian) = compose(left, right)

bloch_hamiltonian_inverse_weight(harmonic::Int) = 1 // harmonic
bloch_liouvillian_inverse_weight(harmonic::Int) = im // harmonic

function floquet_expansion_impl(
  generator::P, gauge::VanVleck{HoriDeprit}, order::Int, provenance::R
) where {P<:PeriodicGenerator,R<:FloquetProvenance}
  return invoke(
    floquet_expansion_impl, Tuple{P,Gauge,Int,R}, generator, gauge, order, provenance
  )
end

function floquet_expansion_impl(
  ::P, ::VanVleck{A}, ::Int, ::R
) where {P<:PeriodicGenerator,A<:ExpansionAlgorithm,R<:FloquetProvenance}
  return throw(
    ArgumentError(
      "no van Vleck expansion implementation for algorithm $(A) and generator $(P)"
    ),
  )
end

function bloch_feshbach_zero_expansion(
  generator::P, gauge::G, order::Int, provenance::R
) where {T,P<:PeriodicGenerator{T},G<:VanVleck{BlochFeshbach},R<:FloquetProvenance}
  kick_components = P[zero(generator) for _ in 1:max(order - 1, 0)]
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

function bloch_phase_connected_log!(
  log_embedding::Vector{Dict{H,T}}, phase::P; simplifier=identity
) where {H,T,P}
  for (n, embedding) in enumerate(log_embedding)
    scale = phase^n
    for harmonic in collect(keys(embedding))
      component = simplifier(scale * embedding[harmonic])::T
      if iszero(component)
        delete!(embedding, harmonic)
      else
        embedding[harmonic] = component
      end
    end
  end
  return log_embedding
end

function evaluate_bloch_connected_van_vleck_phased(
  projection_plan::BlochProjectionPlan{H},
  bloch::BlochProjectionResult{H,T},
  components::AbstractDict{H,T},
  connected_plan::BlochConnectedVanVleckPlan{H},
  log_phase::L;
  product,
  zero_component::T,
  simplifier=identity,
) where {H,T,L}
  canonical_order = length(bloch.effective) - 1
  projection_plan.order == canonical_order + 1 ||
    throw(ArgumentError("Bloch plan/result truncations are inconsistent"))
  length(connected_plan.log.outputs) == canonical_order ||
    throw(ArgumentError("connected-log plan/result truncations are inconsistent"))
  connected_plan.static.zero_harmonic == projection_plan.zero_harmonic ||
    throw(ArgumentError("static-exponential and Bloch zero harmonics are inconsistent"))

  log_embedding = evaluate_bloch_connected_log_plan(
    connected_plan.log, components; product, zero_component, simplifier
  )
  bloch_phase_connected_log!(log_embedding, log_phase; simplifier)

  identity_component = one(first(bloch.effective))
  static_factor = evaluate_bloch_static_exp_plan(
    connected_plan.static,
    log_embedding,
    identity_component,
    zero_component;
    product,
    simplifier,
  )
  counts = BlochVanVleckCounts()
  counts.harmonic_products =
    bloch_connected_log_products(connected_plan.log) + connected_plan.static.product_count
  inverse_static_factor = bloch_vv_static_series_inverse(
    static_factor, canonical_order, product, counts; simplifier
  )
  right_transformed = bloch_vv_static_series_product(
    bloch.effective, static_factor, canonical_order, product, counts; simplifier
  )
  effective = bloch_vv_static_series_product(
    inverse_static_factor, right_transformed, canonical_order, product, counts; simplifier
  )

  return BlochConnectedVanVleckResult(
    static_factor, inverse_static_factor, log_embedding, effective, counts
  )
end

function bloch_feshbach_expansion_impl(
  generator::P,
  gauge::G,
  order::Int,
  provenance::R,
  product::F,
  inverse_weight::W,
  log_phase::L,
  kick_phase::K,
) where {T,P<:PeriodicGenerator{T},G<:VanVleck{BlochFeshbach},R<:FloquetProvenance,F,W,L,K}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  order == 1 && return bloch_feshbach_leading_expansion(generator, gauge, provenance)
  isempty(generator) &&
    return bloch_feshbach_zero_expansion(generator, gauge, order, provenance)

  components = getfield(generator, :components)
  zero_component = getfield(generator, :zero_component)
  projection_plan = compile_bloch_projection_plan(collect(keys(generator)), order)
  bloch = evaluate_bloch_projection_plan(
    projection_plan,
    components;
    product,
    inverse_weight,
    zero_component,
    simplifier=SQA.simplify,
  )
  connected_plan = compile_bloch_connected_van_vleck_plan(projection_plan)
  converted = evaluate_bloch_connected_van_vleck_phased(
    projection_plan,
    bloch,
    components,
    connected_plan,
    log_phase;
    product,
    zero_component,
    simplifier=SQA.simplify,
  )

  kick_components = P[]
  for embedding in converted.log_embedding
    periodic = periodic_generator(embedding, generator.wd, zero_component)
    push!(kick_components, SQA.simplify(kick_phase * periodic)::P)
  end
  effective_components = T[SQA.simplify(component)::T for component in converted.effective]

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

function floquet_expansion_impl(
  generator::PeriodicGenerator{SQA.QAdd},
  gauge::VanVleck{BlochFeshbach},
  order::Int,
  provenance::R,
) where {R<:FloquetProvenance}
  require_hermitian_drive(generator)
  return bloch_feshbach_expansion_impl(
    generator,
    gauge,
    order,
    provenance,
    bloch_hamiltonian_product,
    bloch_hamiltonian_inverse_weight,
    1,
    im,
  )
end

function floquet_expansion_impl(
  generator::PeriodicGenerator{Liouvillian},
  gauge::VanVleck{BlochFeshbach},
  order::Int,
  provenance::R,
) where {R<:FloquetProvenance}
  return bloch_feshbach_expansion_impl(
    generator,
    gauge,
    order,
    provenance,
    bloch_liouvillian_product,
    bloch_liouvillian_inverse_weight,
    im,
    1,
  )
end

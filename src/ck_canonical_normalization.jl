struct CKCanonicalNormalizationCounts
  normalization_products::Int
  inverse_products::Int
  similarity_products::Int
end

ck_canonical_total_products(counts::CKCanonicalNormalizationCounts) =
  counts.normalization_products + counts.inverse_products + counts.similarity_products

struct CKCanonicalNormalizationResult{K}
  static_factor::Vector{K}
  inverse_static_factor::Vector{K}
  log_embedding::Vector{K}
  wave::Vector{K}
  effective::Vector{K}
  counts::CKCanonicalNormalizationCounts
end

ck_endpoint_output_number(key::CKEndpointKey) = count(ck_is_jump, key.vertices)

function ck_endpoint_scale(state::CKEndpointKernel{H,T}, factor) where {H,T}
  result = Dict{CKEndpointKey{H},T}()
  for (key, value) in state.terms
    scaled = convert(T, factor * value)
    ck_endpoint_accumulate!(result, key, scaled, state.zero_component)
  end
  return CKEndpointKernel(result, state.zero_component)
end

function ck_endpoint_project_model(state::CKEndpointKernel{H,T}) where {H,T}
  result = Dict{CKEndpointKey{H},T}()
  for (key, value) in state.terms
    vertex_count = length(key.vertices)
    if iszero(vertex_count)
      ck_endpoint_accumulate!(result, key, value, state.zero_component)
      continue
    end

    whole = CKEndpointConstraint(0, vertex_count)
    whole in key.resolvent_intervals && continue
    model_intervals = copy(key.model_intervals)
    whole in model_intervals || push!(model_intervals, whole)
    projected = CKEndpointKey(
      copy(key.vertices), model_intervals, copy(key.resolvent_intervals)
    )
    ck_endpoint_accumulate!(result, projected, value, state.zero_component)
  end
  return CKEndpointKernel(result, state.zero_component)
end

function evaluate_ck_endpoint_canonical_normalization(
  effective::AbstractVector{K},
  wave::AbstractVector{K},
  order::Int,
  identity_state::K,
  zero_state::K,
) where {K<:CKEndpointKernel}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  length(effective) >= order ||
    throw(ArgumentError("effective series does not contain the requested order"))
  length(wave) >= order - 1 ||
    throw(ArgumentError("wave series does not contain the required canonical orders"))
  identity_state.zero_component == zero_state.zero_component ||
    throw(ArgumentError("identity and zero endpoint states must share one component type"))

  canonical_order = order - 1
  static_factor = K[zero_state for _ in 1:canonical_order]
  inverse_static_factor = K[zero_state for _ in 1:canonical_order]
  normalized_wave = K[zero_state for _ in 1:canonical_order]
  log_embedding = K[zero_state for _ in 1:canonical_order]
  powers = [K[zero_state for _ in 1:canonical_order] for _ in 1:canonical_order]

  normalization_products = 0
  for n in 1:canonical_order
    prefactor = wave[n]
    for j in 1:(n - 1)
      prefactor += ck_endpoint_product(wave[j], static_factor[n - j])
      normalization_products += 1
    end

    nonlinear_log = zero_state
    for power in 2:n
      coefficient = zero_state
      for k in 1:(n - power + 1)
        coefficient += ck_endpoint_product(normalized_wave[k], powers[power - 1][n - k])
        normalization_products += 1
      end
      powers[power][n] = coefficient
      weight = (-1)^(power + 1) // power
      nonlinear_log += ck_endpoint_scale(coefficient, weight)
    end

    candidate = prefactor + nonlinear_log
    static_factor[n] = ck_endpoint_scale(ck_endpoint_project_model(candidate), -1)
    normalized_wave[n] = prefactor + static_factor[n]
    powers[1][n] = normalized_wave[n]
    log_embedding[n] = normalized_wave[n] + nonlinear_log
  end

  inverse_products = 0
  for n in 1:canonical_order
    coefficient = static_factor[n]
    for j in 1:(n - 1)
      coefficient += ck_endpoint_product(static_factor[j], inverse_static_factor[n - j])
      inverse_products += 1
    end
    inverse_static_factor[n] = ck_endpoint_scale(coefficient, -1)
  end

  right_transformed = K[zero_state for _ in 1:order]
  canonical_effective = K[zero_state for _ in 1:order]
  similarity_products = 0
  for n in 1:order
    coefficient = effective[n]
    for j in 1:(n - 1)
      coefficient += ck_endpoint_product(effective[j], static_factor[n - j])
      similarity_products += 1
    end
    right_transformed[n] = coefficient
  end

  for n in 1:order
    coefficient = right_transformed[n]
    for j in 1:(n - 1)
      coefficient += ck_endpoint_product(
        inverse_static_factor[j], right_transformed[n - j]
      )
      similarity_products += 1
    end
    canonical_effective[n] = coefficient
  end

  counts = CKCanonicalNormalizationCounts(
    normalization_products, inverse_products, similarity_products
  )
  return CKCanonicalNormalizationResult(
    static_factor,
    inverse_static_factor,
    log_embedding,
    normalized_wave,
    canonical_effective,
    counts,
  )
end

function evaluate_ck_canonical_normalization(
  effective::AbstractVector{K},
  wave::AbstractVector{K},
  order::Int,
  identity_state::K,
  zero_state::K,
) where {K<:CKPhysicalKernel}
  endpoint_effective = [ck_endpoint_kernel(effective[n]) for n in 1:order]
  endpoint_wave = [ck_endpoint_kernel(wave[n]) for n in 1:(order - 1)]
  identity_endpoint = ck_endpoint_kernel(identity_state)
  zero_endpoint = ck_endpoint_kernel(zero_state)
  return evaluate_ck_endpoint_canonical_normalization(
    endpoint_effective, endpoint_wave, order, identity_endpoint, zero_endpoint
  )
end

function evaluate_ck_period_amplitude(
  effective::AbstractVector{K},
  wave::AbstractVector{K},
  order::Int,
  identity_state::K,
  zero_state::K,
) where {K<:CKEndpointKernel}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  length(effective) >= order ||
    throw(ArgumentError("effective series does not contain the requested order"))
  length(wave) >= order - 1 ||
    throw(ArgumentError("wave series does not contain the required endpoint orders"))

  zero_period = ck_period_zero(zero_state)
  identity_period = ck_period_constant(identity_state)

  endpoint_wave = typeof(identity_period)[zero_period for _ in 0:order]
  endpoint_wave[1] = identity_period
  for n in 1:(order - 1)
    endpoint_wave[n + 1] = ck_period_constant(wave[n])
  end

  endpoint_wave_inverse = ck_unit_series_inverse(
    endpoint_wave, order, identity_period, zero_period
  )

  slow_generator = typeof(identity_period)[zero_period for _ in 0:order]
  for n in 1:order
    slow_generator[n + 1] = ck_period_monomial(effective[n], 1)
  end
  slow_propagator = ck_zero_constant_series_exponential(
    slow_generator, order, identity_period, zero_period
  )

  dressed_left = ck_truncated_series_product(
    endpoint_wave, slow_propagator, order, zero_period
  )
  amplitude = ck_truncated_series_product(
    dressed_left, endpoint_wave_inverse, order, zero_period
  )
  return CKPeriodAmplitudeResult(
    amplitude, endpoint_wave, endpoint_wave_inverse, slow_propagator
  )
end

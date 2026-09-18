struct BlochVanVleckCounts
  factor_products::Int
  log_products::Int
  inverse_products::Int
  similarity_products::Int
  harmonic_products::Int
end

BlochVanVleckCounts() = BlochVanVleckCounts(0, 0, 0, 0, 0)

struct BlochVanVleckResult{H,T}
  static_factor::Vector{T}
  inverse_static_factor::Vector{T}
  normalized_embedding::Vector{Dict{H,T}}
  log_embedding::Vector{Dict{H,T}}
  effective::Vector{T}
  counts::BlochVanVleckCounts
end

function bloch_vv_accumulate!(
  destination::Dict{H,T}, harmonic::H, value::T
) where {H,T}
  destination[harmonic] = haskey(destination, harmonic) ? destination[harmonic] + value : value
  return destination
end

function bloch_vv_simplify_embedding(
  embedding::Dict{H,T}, simplifier
) where {H,T}
  simplified = Dict{H,T}()
  for (harmonic, value) in embedding
    component = simplifier(value)::T
    iszero(component) || (simplified[harmonic] = component)
  end
  return simplified
end

function bloch_vv_add(
  left::Dict{H,T}, right::Dict{H,T}, simplifier
) where {H,T}
  result = copy(left)
  for (harmonic, value) in right
    bloch_vv_accumulate!(result, harmonic, value)
  end
  return bloch_vv_simplify_embedding(result, simplifier)
end

function bloch_vv_scale(
  weight, embedding::Dict{H,T}, simplifier
) where {H,T}
  result = Dict{H,T}()
  for (harmonic, value) in embedding
    component = simplifier(weight * value)::T
    iszero(component) || (result[harmonic] = component)
  end
  return result
end

function bloch_vv_periodic_product(
  left::Dict{H,T},
  right::Dict{H,T},
  product,
  simplifier,
  counts::BlochVanVleckCounts,
) where {H,T}
  result = Dict{H,T}()
  for (left_harmonic, left_value) in left, (right_harmonic, right_value) in right
    counts.harmonic_products += 1
    harmonic = left_harmonic + right_harmonic
    bloch_vv_accumulate!(result, harmonic, product(left_value, right_value))
  end
  return bloch_vv_simplify_embedding(result, simplifier)
end

function bloch_vv_right_static_product(
  periodic::Dict{H,T},
  static::T,
  product,
  simplifier,
  counts::BlochVanVleckCounts,
) where {H,T}
  result = Dict{H,T}()
  for (harmonic, value) in periodic
    counts.harmonic_products += 1
    component = simplifier(product(value, static))::T
    iszero(component) || (result[harmonic] = component)
  end
  return result
end

function bloch_vv_static_series_inverse(
  factor::Vector{T},
  order::Int,
  product,
  counts::BlochVanVleckCounts;
  simplifier=identity,
) where {T}
  inverse = T[first(factor)]
  for n in 1:order
    coefficient = zero(first(factor))
    for k in 1:n
      coefficient += product(factor[k + 1], inverse[n - k + 1])
      counts.inverse_products += 1
    end
    push!(inverse, simplifier(-coefficient)::T)
  end
  return inverse
end

function bloch_vv_static_series_product(
  left::Vector{T},
  right::Vector{T},
  order::Int,
  product,
  counts::BlochVanVleckCounts;
  simplifier=identity,
) where {T}
  result = Vector{T}(undef, order + 1)
  for n in 0:order
    coefficient = zero(first(left))
    for k in 0:n
      coefficient += product(left[k + 1], right[n - k + 1])
      counts.similarity_products += 1
    end
    result[n + 1] = simplifier(coefficient)::T
  end
  return result
end

function bloch_van_vleck_reconstruction(
  plan::BlochProjectionPlan{H},
  bloch::BlochProjectionResult{H,T};
  product,
  zero_component::T,
  simplifier=identity,
) where {H,T}
  order = length(bloch.effective) - 1
  length(bloch.wave) == order ||
    throw(ArgumentError("Bloch wave/effective truncations are inconsistent"))
  plan.order == order + 1 ||
    throw(ArgumentError("Bloch plan/result truncations are inconsistent"))

  identity_component = one(first(bloch.effective))
  static_factor = T[identity_component]
  normalized_embedding = Vector{Dict{H,T}}()
  log_embedding = Vector{Dict{H,T}}()
  powers = [
    [Dict{H,T}() for _ in 1:max(order, 1)] for _ in 1:max(order, 1)
  ]
  counts = BlochVanVleckCounts()

  for n in 1:order
    prefactor = copy(bloch.wave[n])
    for j in 1:(n - 1)
      counts.factor_products += 1
      correction = bloch_vv_right_static_product(
        bloch.wave[j],
        static_factor[n - j + 1],
        product,
        simplifier,
        counts,
      )
      prefactor = bloch_vv_add(prefactor, correction, simplifier)
    end

    nonlinear_log = Dict{H,T}()
    for power in 2:n
      power_coefficient = Dict{H,T}()
      for k in 1:(n - power + 1)
        counts.log_products += 1
        contribution = bloch_vv_periodic_product(
          normalized_embedding[k], powers[power - 1][n - k], product, simplifier, counts
        )
        power_coefficient = bloch_vv_add(power_coefficient, contribution, simplifier)
      end
      powers[power][n] = power_coefficient
      weight = (-1)^(power + 1) * (1 // power)
      nonlinear_log = bloch_vv_add(
        nonlinear_log, bloch_vv_scale(weight, power_coefficient, simplifier), simplifier
      )
    end

    candidate = bloch_vv_add(prefactor, nonlinear_log, simplifier)
    static_n = simplifier(-get(candidate, plan.zero_harmonic, zero_component))::T
    push!(static_factor, static_n)

    normalized_n = copy(prefactor)
    if !iszero(static_n)
      bloch_vv_accumulate!(normalized_n, plan.zero_harmonic, static_n)
      normalized_n = bloch_vv_simplify_embedding(normalized_n, simplifier)
    end
    push!(normalized_embedding, normalized_n)
    powers[1][n] = normalized_n
    push!(log_embedding, bloch_vv_add(normalized_n, nonlinear_log, simplifier))
  end

  inverse_static_factor = bloch_vv_static_series_inverse(
    static_factor, order, product, counts; simplifier
  )
  right_transformed = bloch_vv_static_series_product(
    bloch.effective, static_factor, order, product, counts; simplifier
  )
  effective = bloch_vv_static_series_product(
    inverse_static_factor, right_transformed, order, product, counts; simplifier
  )

  return BlochVanVleckResult(
    static_factor,
    inverse_static_factor,
    normalized_embedding,
    log_embedding,
    effective,
    counts,
  )
end

function bloch_van_vleck_mercator_products(order::Int)
  order >= 0 || throw(ArgumentError("order must be nonnegative"))
  return order < 2 ? 0 : (order + 1) * order * (order - 1) ÷ 6
end

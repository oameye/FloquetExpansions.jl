using FloquetExpansions

include(joinpath(@__DIR__, "bloch_feshbach_reference.jl"))

mutable struct BlochVanVleckReferenceCounts
  factor_products::Int
  log_products::Int
  inverse_products::Int
  similarity_products::Int
  periodic_component_products::Int
end

BlochVanVleckReferenceCounts() = BlochVanVleckReferenceCounts(0, 0, 0, 0, 0)

struct BlochVanVleckReferenceResult{P,T}
  static_factor::Vector{T}
  inverse_static_factor::Vector{T}
  log_embedding::Vector{P}
  effective::Vector{T}
  counts::BlochVanVleckReferenceCounts
end

function counted_periodic_product(
  left::P, right::P, product, counts::BlochVanVleckReferenceCounts
) where {P<:PeriodicGenerator}
  component_counts = BlochReferenceCounts()
  result = reference_periodic_product(left, right, product, component_counts)
  counts.periodic_component_products += component_counts.component_products
  return result
end

function counted_right_static_product(
  periodic::PeriodicGenerator{T}, static::T, product, counts::BlochVanVleckReferenceCounts
) where {T}
  component_counts = BlochReferenceCounts()
  result = reference_right_static_product(periodic, static, product, component_counts)
  counts.periodic_component_products += component_counts.component_products
  return result
end

function reference_static_series_product(
  left::Vector{T},
  right::Vector{T},
  order::Int,
  product,
  counts::BlochVanVleckReferenceCounts;
  simplifier=identity,
) where {T}
  result = Vector{T}(undef, order + 1)
  for n in 0:order
    coefficient = zero(first(left))
    for k in 0:n
      coefficient += product(left[k + 1], right[n - k + 1])
      counts.similarity_products += 1
    end
    result[n + 1] = simplifier(coefficient)
  end
  return result
end

function reference_static_series_inverse(
  factor::Vector{T},
  order::Int,
  product,
  counts::BlochVanVleckReferenceCounts;
  simplifier=identity,
) where {T}
  inverse = T[first(factor)]
  for n in 1:order
    coefficient = zero(first(factor))
    for k in 1:n
      coefficient += product(factor[k + 1], inverse[n - k + 1])
      counts.inverse_products += 1
    end
    push!(inverse, simplifier(-coefficient))
  end
  return inverse
end

function bloch_van_vleck_reference(
  bloch::BlochReferenceResult{P,T}, template::P; product, simplifier=identity
) where {T,P<:PeriodicGenerator{T}}
  order = length(bloch.effective) - 1
  length(bloch.wave) == order ||
    throw(ArgumentError("Bloch wave/effective truncations are inconsistent"))

  identity_component = one(first(bloch.effective))
  static_factor = T[identity_component]
  log_embedding = P[]
  wave_product = P[]
  powers = [[zero(template) for _ in 1:order] for _ in 1:order]
  counts = BlochVanVleckReferenceCounts()

  for n in 1:order
    prefactor = bloch.wave[n]
    for j in 1:(n - 1)
      counts.factor_products += 1
      prefactor += counted_right_static_product(
        bloch.wave[j], static_factor[n - j + 1], product, counts
      )
    end
    prefactor = simplifier(prefactor)

    nonlinear_log = zero(template)
    for power in 2:n
      power_coefficient = zero(template)
      for k in 1:(n - power + 1)
        counts.log_products += 1
        power_coefficient += counted_periodic_product(
          wave_product[k], powers[power - 1][n - k], product, counts
        )
      end
      power_coefficient = simplifier(power_coefficient)
      powers[power][n] = power_coefficient
      weight = (-1)^(power + 1) * (1 // power)
      nonlinear_log += weight * power_coefficient
    end
    nonlinear_log = simplifier(nonlinear_log)

    candidate = simplifier(prefactor + nonlinear_log)
    static_n = simplifier(-time_average(candidate))
    push!(static_factor, static_n)

    wave_n = simplifier(prefactor + reference_static_generator(static_n, template))
    push!(wave_product, wave_n)
    powers[1][n] = wave_n

    generator_n = simplifier(wave_n + nonlinear_log)
    push!(log_embedding, generator_n)
  end

  inverse_static_factor = reference_static_series_inverse(
    static_factor, order, product, counts; simplifier
  )
  right_transformed = reference_static_series_product(
    bloch.effective, static_factor, order, product, counts; simplifier
  )
  effective = reference_static_series_product(
    inverse_static_factor, right_transformed, order, product, counts; simplifier
  )

  return BlochVanVleckReferenceResult(
    static_factor, inverse_static_factor, log_embedding, effective, counts
  )
end

function expected_mercator_products(order::Int)
  order >= 0 || throw(ArgumentError("order must be nonnegative"))
  return order < 2 ? 0 : (order + 1) * order * (order - 1) ÷ 6
end

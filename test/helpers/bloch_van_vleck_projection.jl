using FloquetExpansions

const FE_BVV = FloquetExpansions

mutable struct BlochVanVleckProjectionCounts
  factor_products::Int
  log_products::Int
  inverse_products::Int
  similarity_products::Int
  periodic_component_products::Int
end

BlochVanVleckProjectionCounts() = BlochVanVleckProjectionCounts(0, 0, 0, 0, 0)

struct PeriodicBlochProjectionResult{P,T}
  wave::Vector{P}
  effective::Vector{T}
end

struct BlochVanVleckProjectionResult{P,T}
  static_factor::Vector{T}
  inverse_static_factor::Vector{T}
  normalized_embedding::Vector{P}
  log_embedding::Vector{P}
  effective::Vector{T}
  counts::BlochVanVleckProjectionCounts
end

function projection_periodic_product(
  left::PeriodicGenerator{T},
  right::PeriodicGenerator{T},
  product,
  counts::BlochVanVleckProjectionCounts,
) where {T}
  isequal(left.wd, right.wd) ||
    throw(ArgumentError("periodic generators use different frequencies"))
  out = Dict{Int,T}()
  for left_harmonic in keys(left), right_harmonic in keys(right)
    counts.periodic_component_products += 1
    harmonic = left_harmonic + right_harmonic
    term = product(left[left_harmonic], right[right_harmonic])
    out[harmonic] = haskey(out, harmonic) ? out[harmonic] + term : term
  end
  return PeriodicGenerator(out, left.wd, left.zero_component)
end

function projection_right_static_product(
  periodic::PeriodicGenerator{T},
  static::T,
  product,
  counts::BlochVanVleckProjectionCounts,
) where {T}
  out = Dict{Int,T}()
  for harmonic in keys(periodic)
    counts.periodic_component_products += 1
    out[harmonic] = product(periodic[harmonic], static)
  end
  return PeriodicGenerator(out, periodic.wd, periodic.zero_component)
end

function projection_static_generator(component::T, template::PeriodicGenerator{T}) where {T}
  return PeriodicGenerator(Dict(0 => component), template.wd, template.zero_component)
end

function projection_static_series_inverse(
  factor::Vector{T},
  order::Int,
  product,
  counts::BlochVanVleckProjectionCounts;
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

function projection_static_series_product(
  left::Vector{T},
  right::Vector{T},
  order::Int,
  product,
  counts::BlochVanVleckProjectionCounts;
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

function periodicize_bloch_projection(
  result::FE_BVV.BlochProjectionResult{H,T}, template::PeriodicGenerator{T}
) where {H,T}
  H === Int || throw(ArgumentError("periodic reconstruction currently expects integer harmonics"))
  wave = PeriodicGenerator{T}[
    PeriodicGenerator(Dict(harmonic => value for (harmonic, value) in coefficient), template.wd, template.zero_component)
    for coefficient in result.wave
  ]
  return PeriodicBlochProjectionResult(wave, result.effective)
end

function evaluate_periodic_bloch_projection(
  generator::PeriodicGenerator{T}, order::Int; product, inverse_weight, simplifier=identity
) where {T}
  components = getfield(generator, :components)
  plan = FE_BVV.compile_bloch_projection_plan(collect(keys(generator)), order)
  raw = FE_BVV.evaluate_bloch_projection_plan(
    plan,
    components;
    product,
    inverse_weight,
    zero_component=getfield(generator, :zero_component),
    simplifier,
  )
  return periodicize_bloch_projection(raw, generator)
end

function bloch_van_vleck_projection(
  bloch::PeriodicBlochProjectionResult{P,T}, template::P; product, simplifier=identity
) where {T,P<:PeriodicGenerator{T}}
  order = length(bloch.effective) - 1
  length(bloch.wave) == order ||
    throw(ArgumentError("Bloch wave/effective truncations are inconsistent"))

  identity_component = one(first(bloch.effective))
  static_factor = T[identity_component]
  log_embedding = P[]
  normalized_embedding = P[]
  powers = [[zero(template) for _ in 1:order] for _ in 1:order]
  counts = BlochVanVleckProjectionCounts()

  for n in 1:order
    prefactor = bloch.wave[n]
    for j in 1:(n - 1)
      counts.factor_products += 1
      prefactor += projection_right_static_product(
        bloch.wave[j], static_factor[n - j + 1], product, counts
      )
    end
    prefactor = simplifier(prefactor)

    nonlinear_log = zero(template)
    for power in 2:n
      power_coefficient = zero(template)
      for k in 1:(n - power + 1)
        counts.log_products += 1
        power_coefficient += projection_periodic_product(
          normalized_embedding[k], powers[power - 1][n - k], product, counts
        )
      end
      power_coefficient = simplifier(power_coefficient)
      powers[power][n] = power_coefficient
      nonlinear_log += ((-1)^(power + 1) * (1 // power)) * power_coefficient
    end
    nonlinear_log = simplifier(nonlinear_log)

    candidate = simplifier(prefactor + nonlinear_log)
    static_n = simplifier(-time_average(candidate))
    push!(static_factor, static_n)

    normalized_n = simplifier(prefactor + projection_static_generator(static_n, template))
    push!(normalized_embedding, normalized_n)
    powers[1][n] = normalized_n
    push!(log_embedding, simplifier(normalized_n + nonlinear_log))
  end

  inverse_static_factor = projection_static_series_inverse(
    static_factor, order, product, counts; simplifier
  )
  right_transformed = projection_static_series_product(
    bloch.effective, static_factor, order, product, counts; simplifier
  )
  effective = projection_static_series_product(
    inverse_static_factor, right_transformed, order, product, counts; simplifier
  )

  return BlochVanVleckProjectionResult(
    static_factor,
    inverse_static_factor,
    normalized_embedding,
    log_embedding,
    effective,
    counts,
  )
end

function replay_mercator_log(
  normalized_embedding::Vector{P}, template::P; product, simplifier=identity
) where {P<:PeriodicGenerator}
  order = length(normalized_embedding)
  counts = BlochVanVleckProjectionCounts()
  powers = [[zero(template) for _ in 1:order] for _ in 1:order]
  log_embedding = P[]
  for n in 1:order
    powers[1][n] = normalized_embedding[n]
    generator_n = normalized_embedding[n]
    for power in 2:n
      coefficient = zero(template)
      for k in 1:(n - power + 1)
        counts.log_products += 1
        coefficient += projection_periodic_product(
          normalized_embedding[k], powers[power - 1][n - k], product, counts
        )
      end
      coefficient = simplifier(coefficient)
      powers[power][n] = coefficient
      generator_n += ((-1)^(power + 1) * (1 // power)) * coefficient
    end
    push!(log_embedding, simplifier(generator_n))
  end
  return log_embedding, counts
end

function replay_static_similarity(
  effective::Vector{T}, static_factor::Vector{T}; product, simplifier=identity
) where {T}
  order = length(effective) - 1
  counts = BlochVanVleckProjectionCounts()
  inverse = projection_static_series_inverse(static_factor, order, product, counts; simplifier)
  right = projection_static_series_product(
    effective, static_factor, order, product, counts; simplifier
  )
  transformed = projection_static_series_product(
    inverse, right, order, product, counts; simplifier
  )
  return transformed, counts
end

expected_projection_mercator_products(order::Int) =
  order < 2 ? 0 : (order + 1) * order * (order - 1) ÷ 6

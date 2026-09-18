using FloquetExpansions

const FE_CAYLEY = FloquetExpansions

mutable struct CayleyLogCounts
  cayley_products::Int
  power_products::Int
  harmonic_products::Int
end

CayleyLogCounts() = CayleyLogCounts(0, 0, 0)

function cayley_embedding_product(
  left::Dict{H,T}, right::Dict{H,T}, product, simplifier, counts::CayleyLogCounts
) where {H,T}
  result = Dict{H,T}()
  for (left_harmonic, left_value) in left, (right_harmonic, right_value) in right
    counts.harmonic_products += 1
    harmonic = left_harmonic + right_harmonic
    value = product(left_value, right_value)
    result[harmonic] = haskey(result, harmonic) ? result[harmonic] + value : value
  end
  simplified = Dict{H,T}()
  for (harmonic, value) in result
    component = simplifier(value)::T
    iszero(component) || (simplified[harmonic] = component)
  end
  return simplified
end

function cayley_embedding_add(
  left::Dict{H,T}, right::Dict{H,T}, simplifier
) where {H,T}
  result = copy(left)
  for (harmonic, value) in right
    result[harmonic] = haskey(result, harmonic) ? result[harmonic] + value : value
  end
  simplified = Dict{H,T}()
  for (harmonic, value) in result
    component = simplifier(value)::T
    iszero(component) || (simplified[harmonic] = component)
  end
  return simplified
end

function cayley_embedding_scale(weight, embedding::Dict{H,T}, simplifier) where {H,T}
  result = Dict{H,T}()
  for (harmonic, value) in embedding
    component = simplifier(weight * value)::T
    iszero(component) || (result[harmonic] = component)
  end
  return result
end

function cayley_series_product(
  left::Vector{Dict{H,T}},
  right::Vector{Dict{H,T}},
  product,
  simplifier,
  counts::CayleyLogCounts,
) where {H,T}
  order = length(left)
  length(right) == order || throw(ArgumentError("series truncations are inconsistent"))
  result = [Dict{H,T}() for _ in 1:order]
  for n in 2:order
    coefficient = Dict{H,T}()
    for k in 1:(n - 1)
      isempty(left[k]) && continue
      isempty(right[n - k]) && continue
      counts.power_products += 1
      term = cayley_embedding_product(left[k], right[n - k], product, simplifier, counts)
      coefficient = cayley_embedding_add(coefficient, term, simplifier)
    end
    result[n] = coefficient
  end
  return result
end

function cayley_log_series(
  normalized::Vector{Dict{H,T}}; product, simplifier=identity
) where {H,T}
  order = length(normalized)
  counts = CayleyLogCounts()
  cayley = [Dict{H,T}() for _ in 1:order]

  # Z = (Omega - I)(Omega + I)^(-1), with Omega = I + A.
  # The exact triangular relation (2I + A) Z = A gives
  #   Z_n = 1/2 (A_n - sum_{k=1}^{n-1} A_k Z_{n-k}).
  for n in 1:order
    correction = Dict{H,T}()
    for k in 1:(n - 1)
      isempty(normalized[k]) && continue
      isempty(cayley[n - k]) && continue
      counts.cayley_products += 1
      term = cayley_embedding_product(
        normalized[k], cayley[n - k], product, simplifier, counts
      )
      correction = cayley_embedding_add(correction, term, simplifier)
    end
    candidate = cayley_embedding_add(
      normalized[n], cayley_embedding_scale(-1, correction, simplifier), simplifier
    )
    cayley[n] = cayley_embedding_scale(1 // 2, candidate, simplifier)
  end

  # log(Omega) = 2 atanh(Z) = 2 (Z + Z^3/3 + Z^5/5 + ...).
  logarithm = [cayley_embedding_scale(2, coefficient, simplifier) for coefficient in cayley]
  order < 3 && return logarithm, counts

  cayley_squared = cayley_series_product(cayley, cayley, product, simplifier, counts)
  odd_power = cayley
  for power in 3:2:order
    odd_power = cayley_series_product(
      odd_power, cayley_squared, product, simplifier, counts
    )
    weight = 2 // power
    for n in power:order
      isempty(odd_power[n]) && continue
      logarithm[n] = cayley_embedding_add(
        logarithm[n], cayley_embedding_scale(weight, odd_power[n], simplifier), simplifier
      )
    end
  end
  return logarithm, counts
end

function mercator_log_series(
  normalized::Vector{Dict{H,T}}; product, simplifier=identity
) where {H,T}
  order = length(normalized)
  counts = CayleyLogCounts()
  powers = [normalized]
  logarithm = [copy(coefficient) for coefficient in normalized]
  previous = normalized
  for power in 2:order
    current = cayley_series_product(previous, normalized, product, simplifier, counts)
    push!(powers, current)
    weight = (-1)^(power + 1) * (1 // power)
    for n in power:order
      isempty(current[n]) && continue
      logarithm[n] = cayley_embedding_add(
        logarithm[n], cayley_embedding_scale(weight, current[n], simplifier), simplifier
      )
    end
    previous = current
  end
  return logarithm, counts
end

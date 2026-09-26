mutable struct BlochVanVleckCounts
  factor_products::Int
  log_products::Int
  inverse_products::Int
  similarity_products::Int
  harmonic_products::Int
end

BlochVanVleckCounts() = BlochVanVleckCounts(0, 0, 0, 0, 0)

function bloch_vv_static_series_inverse(
  factor::Vector{T}, order::Int, product, counts::BlochVanVleckCounts; simplifier=identity
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

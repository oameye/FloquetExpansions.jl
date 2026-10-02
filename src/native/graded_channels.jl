struct GradedChannel{T}
  onset::Int
  weight::T
  coefficients::Vector{Vector{T}}
end

function GradedChannel{T}(onset::Int, coefficients::Vector{Vector{T}}) where {T}
  return GradedChannel{T}(onset, one(T), coefficients)
end

function active_weights(
  channels::AbstractVector{GradedChannel{T}}, indices::AbstractVector{Int}
) where {T}
  return T[channels[index].weight for index in indices]
end

function active_channels(
  channels::AbstractVector{GradedChannel{T}}, order::Int, dimension::Int
) where {T}
  indices = Int[]
  for (index, channel) in pairs(channels)
    channel.onset < order && push!(indices, index)
  end
  active = Matrix{T}(undef, dimension, length(indices))
  for (column, index) in pairs(indices)
    active[:, column] = channels[index].coefficients[1]
  end
  return indices, active
end

function hermitian_part(X::AbstractMatrix)
  S = typeof((zero(eltype(X)) + zero(eltype(X))) / 2)
  result = Matrix{S}(undef, size(X))
  for j in axes(X, 2), i in axes(X, 1)
    result[i, j] = (X[i, j] + conj(X[j, i])) / 2
  end
  return result
end

# result += weight * u * v'
function add_outer!(
  result::AbstractMatrix{T}, weight::T, u::AbstractVector{T}, v::AbstractVector{T}
) where {T}
  for j in eachindex(v)
    conjugate = conj(v[j])
    for i in eachindex(u)
      result[i, j] += (weight * u[i]) * conjugate
    end
  end
  return result
end

function known_gram(
  channels::AbstractVector{GradedChannel{T}}, order::Int, dimension::Int
) where {T}
  result = exact_zeros(T, dimension, dimension)
  for channel in channels
    q = order - channel.onset
    q >= 2 || continue
    length(channel.coefficients) >= q ||
      throw(ArgumentError("graded channel is missing a lower-order amplitude"))
    for k in 1:(q - 1)
      add_outer!(
        result, channel.weight, channel.coefficients[k + 1], channel.coefficients[q - k + 1]
      )
    end
  end
  return hermitian_part(result)
end

function gram_coefficient(
  channels::AbstractVector{GradedChannel{T}}, order::Int, dimension::Int
) where {T}
  result = exact_zeros(T, dimension, dimension)
  for channel in channels
    q = order - channel.onset
    stored = length(channel.coefficients)
    for k in max(0, q + 1 - stored):min(q, stored - 1)
      add_outer!(
        result, channel.weight, channel.coefficients[k + 1], channel.coefficients[q - k + 1]
      )
    end
  end
  return hermitian_part(result)
end

function store_corrections!(
  channels::AbstractVector{GradedChannel{T}},
  indices::AbstractVector{Int},
  correction::AbstractMatrix,
  order::Int,
) where {T}
  size(correction, 2) == length(indices) ||
    throw(DimensionMismatch("one correction column is required per active channel"))
  for (column, index) in pairs(indices)
    channel = channels[index]
    length(channel.coefficients) == order - channel.onset ||
      throw(ArgumentError("graded channel coefficient sequence is not prefix complete"))
    push!(channel.coefficients, Vector{T}(correction[:, column]))
  end
  return channels
end

function store_births!(
  channels::AbstractVector{GradedChannel{T}},
  newborn::AbstractMatrix,
  weights::AbstractVector,
  order::Int,
) where {T}
  length(weights) == size(newborn, 2) ||
    throw(DimensionMismatch("one rate weight is required per newborn channel"))
  for column in axes(newborn, 2)
    push!(
      channels, GradedChannel{T}(order, T(weights[column]), [Vector{T}(newborn[:, column])])
    )
  end
  return channels
end

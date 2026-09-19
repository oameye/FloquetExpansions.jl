struct StinespringHistory{K,P,T}
  output_number::Int
  physical_key::K
  history_path::P
  coefficients::Vector{T}
end

struct StinespringRow{K,T}
  output_number::Int
  physical_key::K
  coefficients::Vector{T}
end

function stinespring_zero_like(value::AbstractMatrix)
  return zero(value)
end

function stinespring_coefficient(
  coefficients::Vector{T}, n::Int, zero_component::T
) where {T}
  return n + 1 <= length(coefficients) ? coefficients[n + 1] : zero_component
end

function coalesce_stinespring_histories(histories::Vector{H}) where {H<:StinespringHistory}
  isempty(histories) && return StinespringRow[]
  zero_component = zero(first(first(histories).coefficients))
  grouped = Dict{Tuple{Int,Any},Vector{typeof(zero_component)}}()

  for history in histories
    key = (history.output_number, history.physical_key)
    coefficients = get!(grouped, key) do
      return typeof(zero_component)[]
    end
    if length(coefficients) < length(history.coefficients)
      append!(
        coefficients,
        [zero_component for _ in 1:(length(history.coefficients) - length(coefficients))],
      )
    end
    for n in eachindex(history.coefficients)
      coefficients[n] += history.coefficients[n]
    end
  end

  rows = StinespringRow[]
  for ((output_number, physical_key), coefficients) in grouped
    push!(rows, StinespringRow(output_number, physical_key, coefficients))
  end
  return rows
end

function stinespring_metric_coefficients(rows, order::Int, identity_component)
  order >= 0 || throw(ArgumentError("order must be nonnegative"))
  zero_component = zero(identity_component)
  metric = [zero_component for _ in 0:order]

  for row in rows
    q = row.output_number
    q < 0 && throw(ArgumentError("output number must be nonnegative"))
    q > order && continue
    for r in q:order
      amplitude_order = r - q
      for left_order in 0:amplitude_order
        right_order = amplitude_order - left_order
        left = stinespring_coefficient(row.coefficients, left_order, zero_component)
        right = stinespring_coefficient(row.coefficients, right_order, zero_component)
        metric[r + 1] += adjoint(left) * right
      end
    end
  end
  return metric
end

function stinespring_series_product(left, right, order::Int, zero_component)
  result = [zero_component for _ in 0:order]
  for n in 0:order
    for k in 0:n
      k + 1 <= length(left) || continue
      n - k + 1 <= length(right) || continue
      result[n + 1] += left[k + 1] * right[n - k + 1]
    end
  end
  return result
end

function stinespring_inverse_series(series, order::Int, identity_component)
  series[1] == identity_component ||
    throw(ArgumentError("reference inverse expects unit leading coefficient"))
  zero_component = zero(identity_component)
  result = [zero_component for _ in 0:order]
  result[1] = identity_component
  for n in 1:order
    coefficient = zero_component
    for k in 1:n
      k + 1 <= length(series) || continue
      coefficient += series[k + 1] * result[n - k + 1]
    end
    result[n + 1] = -coefficient
  end
  return result
end

function stinespring_inverse_sqrt_series(metric, order::Int, identity_component)
  metric[1] == identity_component ||
    throw(ArgumentError("reference inverse square root expects unit leading metric"))
  zero_component = zero(identity_component)
  result = [zero_component for _ in 0:order]
  result[1] = identity_component

  for n in 1:order
    lower = zero_component
    for left_order in 0:(n - 1)
      left_order + 1 <= length(result) || continue
      for metric_order in 0:(n - left_order)
        metric_order + 1 <= length(metric) || continue
        right_order = n - left_order - metric_order
        right_order < n || continue
        right_order + 1 <= length(result) || continue
        lower += result[left_order + 1] * metric[metric_order + 1] * result[right_order + 1]
      end
    end
    result[n + 1] = -(1 // 2) * lower
  end
  return result
end

function stinespring_normalize_rows_series(
  rows, inverse_sqrt, order::Int, identity_component
)
  zero_component = zero(identity_component)
  normalized = StinespringRow[]
  for row in rows
    coefficients = [zero_component for _ in 0:order]
    for n in 0:order
      for k in 0:n
        k + 1 <= length(row.coefficients) || continue
        n - k + 1 <= length(inverse_sqrt) || continue
        coefficients[n + 1] += row.coefficients[k + 1] * inverse_sqrt[n - k + 1]
      end
    end
    push!(normalized, StinespringRow(row.output_number, row.physical_key, coefficients))
  end
  return normalized
end

function stinespring_triangle_contains(
  generator_order::Int, output_number::Int, amplitude_order::Int
)
  generator_order >= 0 || throw(ArgumentError("generator order must be nonnegative"))
  output_number >= 0 || throw(ArgumentError("output number must be nonnegative"))
  amplitude_order >= 0 || throw(ArgumentError("amplitude order must be nonnegative"))
  return output_number + amplitude_order <= generator_order + 1
end

function stinespring_cayley_rows(generator_order::Int, unitary)
  identity_component = Matrix(unitary' * unitary)
  zero_component = zero(identity_component)

  no_jump = [zero_component for _ in 0:(generator_order + 1)]
  no_jump[1] = identity_component
  for n in 1:(generator_order + 1)
    no_jump[n + 1] = (2 * (-1)^n) * identity_component
  end

  one_jump = [zero_component for _ in 0:generator_order]
  for n in 0:generator_order
    one_jump[n + 1] = (2 * (-1)^n) * unitary
  end

  return [StinespringRow(0, :vacuum, no_jump), StinespringRow(1, :one_jump, one_jump)]
end

function stinespring_evaluate_row(row::StinespringRow, T)
  value = zero(first(row.coefficients))
  for (index, coefficient) in enumerate(row.coefficients)
    value += T^(index - 1) * coefficient
  end
  return value
end

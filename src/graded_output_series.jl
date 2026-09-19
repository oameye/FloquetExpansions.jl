struct GradedOutputSeries{K,T}
  coefficients::Vector{Dict{K,T}}
  zero_component::T
end

function graded_output_series(::Type{K}, zero_component::T, cutoff::Int) where {K,T}
  cutoff >= 0 || throw(ArgumentError("cutoff must be nonnegative"))
  return GradedOutputSeries{K,T}(
    [Dict{K,T}() for _ in 0:cutoff], zero_component
  )
end

graded_output_cutoff(series::GradedOutputSeries) = length(series.coefficients) - 1

function graded_output_accumulate!(
  series::GradedOutputSeries{K,T}, degree::Int, key::K, value::T
) where {K,T}
  0 <= degree <= graded_output_cutoff(series) ||
    throw(ArgumentError("degree $(degree) is outside the retained graded series"))
  destination = series.coefficients[degree + 1]
  updated = get(destination, key, series.zero_component) + value
  if updated == series.zero_component
    haskey(destination, key) && delete!(destination, key)
  else
    destination[key] = updated
  end
  return series
end

function graded_output_require_compatible(
  left::GradedOutputSeries{K,T}, right::GradedOutputSeries{K,T}
) where {K,T}
  graded_output_cutoff(left) == graded_output_cutoff(right) ||
    throw(ArgumentError("graded-series cutoffs must agree"))
  left.zero_component == right.zero_component ||
    throw(ArgumentError("graded-series zero components must agree"))
  return nothing
end

function Base.:+(
  left::GradedOutputSeries{K,T}, right::GradedOutputSeries{K,T}
) where {K,T}
  graded_output_require_compatible(left, right)
  result = graded_output_series(K, left.zero_component, graded_output_cutoff(left))
  for degree in 0:graded_output_cutoff(left)
    for (key, value) in left.coefficients[degree + 1]
      graded_output_accumulate!(result, degree, key, value)
    end
    for (key, value) in right.coefficients[degree + 1]
      graded_output_accumulate!(result, degree, key, value)
    end
  end
  return result
end

function Base.:-(
  left::GradedOutputSeries{K,T}, right::GradedOutputSeries{K,T}
) where {K,T}
  return left + (-one(Int)) * right
end

function Base.:*(weight::Number, series::GradedOutputSeries{K,T}) where {K,T}
  result = graded_output_series(K, series.zero_component, graded_output_cutoff(series))
  for degree in 0:graded_output_cutoff(series)
    for (key, value) in series.coefficients[degree + 1]
      scaled = (weight * value)::T
      scaled == series.zero_component ||
        graded_output_accumulate!(result, degree, key, scaled)
    end
  end
  return result
end

Base.:*(series::GradedOutputSeries, weight::Number) = weight * series

struct OutputSectorComposition{K,W}
  key::K
  weight::W
end

OutputSectorComposition(key::K) where {K} = OutputSectorComposition(key, one(Int))

struct GradedOutputProduct{F,G}
  compose_output::F
  system_product::G
end

function (algebra::GradedOutputProduct)(
  left::GradedOutputSeries{K,T}, right::GradedOutputSeries{K,T}
) where {K,T}
  graded_output_require_compatible(left, right)
  cutoff = graded_output_cutoff(left)
  result = graded_output_series(K, left.zero_component, cutoff)

  for left_degree in 0:cutoff
    left_terms = left.coefficients[left_degree + 1]
    isempty(left_terms) && continue
    for right_degree in 0:(cutoff - left_degree)
      right_terms = right.coefficients[right_degree + 1]
      isempty(right_terms) && continue
      degree = left_degree + right_degree
      for (left_key, left_value) in left_terms,
        (right_key, right_value) in right_terms

        composition = algebra.compose_output(left_key, right_key)
        output_key = composition.key::K
        system_value = algebra.system_product(left_value, right_value)
        value = (composition.weight * system_value)::T
        graded_output_accumulate!(result, degree, output_key, value)
      end
    end
  end
  return result
end

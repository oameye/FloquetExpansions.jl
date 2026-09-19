struct CKOutputLeg{C,S}
  channel::C
  sideband::S
end

function ck_output_key(path::Tuple)
  isempty(path) && return ()
  ordered = sort(collect(path); by=leg -> (string(leg.channel), leg.sideband))
  return Tuple(ordered)
end

ck_compose_output(left::Tuple, right::Tuple) = ck_output_key((left..., right...))

struct CKPhysicalVertex{C,H,S,T}
  order::Int
  channel::C
  system_harmonic::H
  output_sideband::S
  mismatch::H
  grade::Int
  value::T
end

function ck_jump_vertex(
  order::Int, channel::C, system_harmonic::H, output_sideband::H, value::T
) where {C,H,T}
  return CKPhysicalVertex(
    order,
    channel,
    system_harmonic,
    output_sideband,
    system_harmonic - output_sideband,
    1,
    value,
  )
end

function ck_drift_vertex(
  order::Int, channel::C, system_harmonic::H, output_sideband::H, value::T
) where {C,H,T}
  return CKPhysicalVertex(
    order, channel, system_harmonic, output_sideband, system_harmonic, 0, value
  )
end

struct CKDeltaOutputSeries{T}
  coefficients::Vector{Dict{Tuple,T}}
  zero_component::T
end

function ck_zero_series(zero_component::T, cutoff::Int) where {T}
  cutoff >= 0 || throw(ArgumentError("cutoff must be nonnegative"))
  return CKDeltaOutputSeries([Dict{Tuple,T}() for _ in 0:cutoff], zero_component)
end

ck_series_cutoff(series::CKDeltaOutputSeries) = length(series.coefficients) - 1

function ck_accumulate!(destination::Dict{Tuple,T}, key::Tuple, value::T, zero_component::T) where {T}
  updated = get(destination, key, zero_component) + value
  if updated == zero_component
    haskey(destination, key) && delete!(destination, key)
  else
    destination[key] = updated
  end
  return destination
end

function Base.:+(left::CKDeltaOutputSeries{T}, right::CKDeltaOutputSeries{T}) where {T}
  cutoff = ck_series_cutoff(left)
  cutoff == ck_series_cutoff(right) || throw(ArgumentError("series cutoffs must agree"))
  result = ck_zero_series(left.zero_component, cutoff)
  for degree in 0:cutoff
    destination = result.coefficients[degree + 1]
    for (key, value) in left.coefficients[degree + 1]
      destination[key] = value
    end
    for (key, value) in right.coefficients[degree + 1]
      ck_accumulate!(destination, key, value, left.zero_component)
    end
  end
  return result
end

function Base.:-(left::CKDeltaOutputSeries{T}, right::CKDeltaOutputSeries{T}) where {T}
  return left + (-one(Int)) * right
end

function Base.:*(weight::Number, series::CKDeltaOutputSeries{T}) where {T}
  cutoff = ck_series_cutoff(series)
  result = ck_zero_series(series.zero_component, cutoff)
  for degree in 0:cutoff
    destination = result.coefficients[degree + 1]
    for (key, value) in series.coefficients[degree + 1]
      scaled = weight * value
      scaled == series.zero_component || (destination[key] = scaled)
    end
  end
  return result
end

function ck_series_product(
  left::CKDeltaOutputSeries{T},
  right::CKDeltaOutputSeries{T},
  products::Base.RefValue{Int},
) where {T}
  cutoff = ck_series_cutoff(left)
  cutoff == ck_series_cutoff(right) || throw(ArgumentError("series cutoffs must agree"))
  result = ck_zero_series(left.zero_component, cutoff)
  for left_degree in 0:cutoff
    left_terms = left.coefficients[left_degree + 1]
    isempty(left_terms) && continue
    for right_degree in 0:(cutoff - left_degree)
      right_terms = right.coefficients[right_degree + 1]
      isempty(right_terms) && continue
      destination = result.coefficients[left_degree + right_degree + 1]
      for (left_key, left_value) in left_terms, (right_key, right_value) in right_terms
        products[] += 1
        output_key = ck_compose_output(left_key, right_key)
        ck_accumulate!(
          destination, output_key, left_value * right_value, left.zero_component
        )
      end
    end
  end
  return result
end

function ck_series_components(vertices_by_order, zero_component::T, cutoff::Int) where {T}
  components = Dict{Int,CKDeltaOutputSeries{T}}()
  for vertices in vertices_by_order, vertex in vertices
    series = get!(components, vertex.mismatch) do
      return ck_zero_series(zero_component, cutoff)
    end
    key = if iszero(vertex.grade)
      ()
    else
      (CKOutputLeg(vertex.channel, vertex.output_sideband),)
    end
    ck_accumulate!(
      series.coefficients[vertex.order + 1], key, vertex.value, zero_component
    )
  end
  return components
end

struct CKDirectState{H}
  harmonic::H
  output::Tuple
end

function ck_direct_tiered_recurrence(
  vertices_by_order,
  order::Int,
  zero_component::T,
  identity_component::T;
  inverse_weight,
) where {T}
  wave = [Dict{CKDirectState{Int},T}() for _ in 1:order]
  effective = [Dict{CKDirectState{Int},T}() for _ in 1:order]
  products = Ref(0)
  identity_state = CKDirectState(0, ())

  for n in 1:order
    residual = Dict{CKDirectState{Int},T}()

    for vertex_order in 1:min(n, length(vertices_by_order))
      previous_order = n - vertex_order
      for vertex in vertices_by_order[vertex_order]
        vertex_key = if iszero(vertex.grade)
          ()
        else
          (CKOutputLeg(vertex.channel, vertex.output_sideband),)
        end
        if iszero(previous_order)
          output = CKDirectState(vertex.mismatch, ck_compose_output(vertex_key, identity_state.output))
          products[] += 1
          value = vertex.value * identity_component
          residual[output] = get(residual, output, zero_component) + value
        else
          for (previous, previous_value) in wave[previous_order]
            output = CKDirectState(
              vertex.mismatch + previous.harmonic,
              ck_compose_output(vertex_key, previous.output),
            )
            products[] += 1
            value = vertex.value * previous_value
            residual[output] = get(residual, output, zero_component) + value
          end
        end
      end
    end

    for wave_order in 1:(n - 1)
      effective_order = n - wave_order
      for (wave_state, wave_value) in wave[wave_order],
        (effective_state, effective_value) in effective[effective_order]

        output = CKDirectState(
          wave_state.harmonic + effective_state.harmonic,
          ck_compose_output(wave_state.output, effective_state.output),
        )
        products[] += 1
        value = wave_value * effective_value
        residual[output] = get(residual, output, zero_component) - value
      end
    end

    for (state, value) in residual
      if iszero(state.harmonic)
        effective[n][state] = value
      else
        wave[n][state] = inverse_weight(state.harmonic) * value
      end
    end
  end

  return (; wave, effective, products=products[])
end

function ck_materialize_bloch_wave(
  bloch, degree::Int, zero_component::T
) where {T}
  result = Dict{CKDirectState{Int},T}()
  for insertion_order in bloch.wave
    for (harmonic, series) in insertion_order
      for (output, value) in series.coefficients[degree + 1]
        state = CKDirectState(harmonic, output)
        result[state] = get(result, state, zero_component) + value
      end
    end
  end
  return Dict(state => value for (state, value) in result if value != zero_component)
end

function ck_materialize_bloch_effective(
  bloch, degree::Int, zero_component::T
) where {T}
  result = Dict{CKDirectState{Int},T}()
  for series in bloch.effective
    for (output, value) in series.coefficients[degree + 1]
      state = CKDirectState(0, output)
      result[state] = get(result, state, zero_component) + value
    end
  end
  return Dict(state => value for (state, value) in result if value != zero_component)
end

function ck_trivial_series_components(
  components::AbstractDict{H,T}, zero_component::T, cutoff::Int
) where {H,T}
  result = Dict{H,CKDeltaOutputSeries{T}}()
  for (harmonic, value) in components
    series = ck_zero_series(zero_component, cutoff)
    series.coefficients[2][()] = value
    result[harmonic] = series
  end
  return result
end

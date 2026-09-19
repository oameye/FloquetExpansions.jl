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

function ck_zero_series(zero_component::T, cutoff::Int) where {T}
  return FloquetExpansions.graded_output_series(Tuple, zero_component, cutoff)
end

function ck_counted_series_product(products::Base.RefValue{Int})
  compose_output =
    (left, right) ->
      FloquetExpansions.OutputSectorComposition(ck_compose_output(left, right))
  system_product = function (left, right)
    products[] += 1
    return left * right
  end
  return FloquetExpansions.GradedOutputProduct(compose_output, system_product)
end

function ck_series_components(vertices_by_order, zero_component::T, cutoff::Int) where {T}
  S = FloquetExpansions.GradedOutputSeries{Tuple,T}
  components = Dict{Int,S}()
  for vertices in vertices_by_order, vertex in vertices
    series = get!(components, vertex.mismatch) do
      return ck_zero_series(zero_component, cutoff)
    end
    key = if iszero(vertex.grade)
      ()
    else
      (CKOutputLeg(vertex.channel, vertex.output_sideband),)
    end
    FloquetExpansions.graded_output_accumulate!(series, vertex.order, key, vertex.value)
  end
  return components
end

struct CKDirectState{H}
  harmonic::H
  output::Tuple
end

function ck_direct_tiered_recurrence(
  vertices_by_order, order::Int, zero_component::T, identity_component::T; inverse_weight
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
          output = CKDirectState(
            vertex.mismatch, ck_compose_output(vertex_key, identity_state.output)
          )
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

function ck_materialize_bloch_wave(bloch, degree::Int, zero_component::T) where {T}
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

function ck_materialize_bloch_effective(bloch, degree::Int, zero_component::T) where {T}
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
  S = FloquetExpansions.GradedOutputSeries{Tuple,T}
  result = Dict{H,S}()
  for (harmonic, value) in components
    series = ck_zero_series(zero_component, cutoff)
    FloquetExpansions.graded_output_accumulate!(series, 1, (), value)
    result[harmonic] = series
  end
  return result
end

if !isdefined(@__MODULE__, :PhysicalOpenLegVertex)
  include(joinpath(@__DIR__, "open_leg_physical_provenance.jl"))
end

struct PersistentOpenLegPathNode
  leg::Int
  tail::Int
end

mutable struct PersistentOpenLegPathRegistry{L}
  legs::Vector{L}
  leg_ids::Dict{L,Int}
  nodes::Vector{PersistentOpenLegPathNode}
  node_ids::Dict{Tuple{Int,Int},Int}
end

function PersistentOpenLegPathRegistry(::Type{L}) where {L}
  return PersistentOpenLegPathRegistry{L}(
    L[], Dict{L,Int}(), [PersistentOpenLegPathNode(0, 0)], Dict{Tuple{Int,Int},Int}()
  )
end

function persistent_openleg_cons!(
  registry::PersistentOpenLegPathRegistry{L}, leg::L, tail::Int
) where {L}
  leg_id = get!(registry.leg_ids, leg) do
    push!(registry.legs, leg)
    return length(registry.legs)
  end
  node_key = (leg_id, tail)
  return get!(registry.node_ids, node_key) do
    push!(registry.nodes, PersistentOpenLegPathNode(leg_id, tail))
    return length(registry.nodes)
  end
end

function persistent_openleg_path!(
  registry::PersistentOpenLegPathRegistry{L}, path::Tuple
) where {L}
  state = 1
  for leg in Iterators.reverse(path)
    state = persistent_openleg_cons!(registry, leg, state)
  end
  return state
end

function persistent_openleg_materialize(
  registry::PersistentOpenLegPathRegistry{L}, state::Int
) where {L}
  legs = L[]
  while state != 1
    node = registry.nodes[state]
    push!(legs, registry.legs[node.leg])
    state = node.tail
  end
  return Tuple(legs)
end

function persistent_openleg_registry(paths::Vector{Tuple})
  for path in paths
    isempty(path) && continue
    return persistent_openleg_registry(paths, typeof(first(path)))
  end
  return throw(ArgumentError("at least one nonempty path is required"))
end

function persistent_openleg_registry(paths::Vector{Tuple}, ::Type{L}) where {L}
  registry = PersistentOpenLegPathRegistry(L)
  ids = Int[]
  for path in paths
    push!(ids, persistent_openleg_path!(registry, path))
  end
  return (; registry, ids)
end

mutable struct PhysicalOpenLegOutputRegistry
  keys::Vector{Tuple}
  ids::Dict{Tuple,Int}
  composition_ids::Dict{Tuple{Int,Int},Int}
end

function PhysicalOpenLegOutputRegistry()
  return PhysicalOpenLegOutputRegistry(
    Tuple[()], Dict{Tuple,Int}(() => 1), Dict{Tuple{Int,Int},Int}()
  )
end

function physical_openleg_intern_output!(
  registry::PhysicalOpenLegOutputRegistry, path::Tuple
)
  key = physical_openleg_output_key(path)
  return get!(registry.ids, key) do
    push!(registry.keys, key)
    return length(registry.keys)
  end
end

function physical_openleg_compose_output!(
  registry::PhysicalOpenLegOutputRegistry, left::Int, right::Int
)
  pair = left <= right ? (left, right) : (right, left)
  return get!(registry.composition_ids, pair) do
    return physical_openleg_intern_output!(
      registry, (registry.keys[left]..., registry.keys[right]...)
    )
  end
end

struct PhysicalOpenLegRepresentationCounts
  flat_paths::Int
  flat_leg_cells::Int
  persistent_nodes::Int
  persistent_leg_atoms::Int
  physical_output_classes::Int
  physical_output_key_cells::Int
end

function physical_openleg_representation_counts(plan::OpenLegEvaluationPlan)
  persistent = persistent_openleg_registry(plan.input.paths)
  for (path, state) in zip(plan.input.paths, persistent.ids)
    persistent_openleg_materialize(persistent.registry, state) == path ||
      error("persistent path round-trip failed")
  end

  outputs = PhysicalOpenLegOutputRegistry()
  for path in plan.input.paths
    physical_openleg_intern_output!(outputs, path)
  end

  return PhysicalOpenLegRepresentationCounts(
    count(path -> !isempty(path), plan.input.paths),
    sum(length(path) for path in plan.input.paths),
    length(persistent.registry.nodes) - 1,
    length(persistent.registry.legs),
    length(outputs.keys) - 1,
    sum(length(key) for key in outputs.keys),
  )
end

struct PhysicalOutputState{H}
  harmonic::H
  output::Int
end

struct PhysicalOutputVertex{H}
  harmonic::H
  output::Int
  component::Int
end

struct PhysicalOutputRecurrenceCounts
  residual_states::Int
  wave_states::Int
  effective_states::Int
  generator_products::Int
  fold_products::Int
  physical_output_classes::Int
  output_compositions::Int
end

function evaluate_physical_output_recurrence(
  vertices_by_order::Vector{Vector{V}},
  order::Int,
  zero_component::T;
  product,
  inverse_weight,
  identity_component::T,
) where {C,T,V<:PhysicalOpenLegVertex{C,Int,Int,T}}
  order >= 1 || throw(ArgumentError("order must be >= 1"))

  outputs = PhysicalOpenLegOutputRegistry()
  values = T[]
  compiled_by_order = Vector{Vector{PhysicalOutputVertex{Int}}}()

  for (vertex_order, vertices) in enumerate(vertices_by_order)
    compiled = PhysicalOutputVertex{Int}[]
    for vertex in vertices
      vertex.order == vertex_order ||
        throw(ArgumentError("vertex order does not match its tier"))
      vertex.grade in (0, 1) ||
        throw(ArgumentError("bare physical reference vertices must have grade zero or one"))
      path = if iszero(vertex.grade)
        ()
      else
        (OpenLegOutputLeg(vertex.channel, vertex.output_sideband),)
      end
      output = physical_openleg_intern_output!(outputs, path)
      push!(values, vertex.value)
      push!(compiled, PhysicalOutputVertex(vertex.mismatch, output, length(values)))
    end
    push!(compiled_by_order, compiled)
  end

  identity_state = PhysicalOutputState(0, 1)
  wave = [Dict{PhysicalOutputState{Int},T}() for _ in 1:order]
  effective = [Dict{PhysicalOutputState{Int},T}() for _ in 1:order]

  residual_state_count = 0
  generator_product_count = 0
  fold_product_count = 0

  for n in 1:order
    residual = Dict{PhysicalOutputState{Int},T}()

    for vertex_order in 1:min(n, length(compiled_by_order))
      previous_order = n - vertex_order
      for vertex in compiled_by_order[vertex_order]
        if iszero(previous_order)
          output = PhysicalOutputState(
            vertex.harmonic,
            physical_openleg_compose_output!(outputs, vertex.output, identity_state.output),
          )
          term = product(values[vertex.component], identity_component)
          residual[output] = get(residual, output, zero_component) + term
          generator_product_count += 1
        else
          for (previous, previous_value) in wave[previous_order]
            output = PhysicalOutputState(
              vertex.harmonic + previous.harmonic,
              physical_openleg_compose_output!(outputs, vertex.output, previous.output),
            )
            term = product(values[vertex.component], previous_value)
            residual[output] = get(residual, output, zero_component) + term
            generator_product_count += 1
          end
        end
      end
    end

    for wave_order in 1:(n - 1)
      effective_order = n - wave_order
      for (wave_state, wave_value) in wave[wave_order],
        (effective_state, effective_value) in effective[effective_order]

        output = PhysicalOutputState(
          wave_state.harmonic + effective_state.harmonic,
          physical_openleg_compose_output!(
            outputs, wave_state.output, effective_state.output
          ),
        )
        term = product(wave_value, effective_value)
        residual[output] = get(residual, output, zero_component) - term
        fold_product_count += 1
      end
    end

    residual_state_count += length(residual)
    for (state, value) in residual
      if iszero(state.harmonic)
        effective[n][state] = value
      else
        wave[n][state] = inverse_weight(state.harmonic) * value
      end
    end
  end

  counts = PhysicalOutputRecurrenceCounts(
    residual_state_count,
    sum(length, wave),
    sum(length, effective),
    generator_product_count,
    fold_product_count,
    length(outputs.keys) - 1,
    length(outputs.composition_ids),
  )
  return (; wave, effective, outputs, counts)
end

function coalesce_flat_physical_outputs(
  components::Dict{OpenLegPlanState{H},T}, plan::OpenLegEvaluationPlan{H,T}
) where {H,T}
  coalesced = Dict{Tuple,T}()
  for (state, value) in components
    output_key = physical_openleg_output_key(plan.input.paths[state.path])
    state.grade == length(output_key) ||
      error("physical output grade is not determined by its key")
    key = (state.harmonic, output_key)
    coalesced[key] = get(coalesced, key, plan.input.zero_component) + value
  end
  return coalesced
end

function materialize_compact_physical_outputs(
  components::Dict{PhysicalOutputState{H},T},
  outputs::PhysicalOpenLegOutputRegistry,
  zero_component::T,
) where {H,T}
  materialized = Dict{Tuple,T}()
  for (state, value) in components
    key = (state.harmonic, outputs.keys[state.output])
    materialized[key] = get(materialized, key, zero_component) + value
  end
  return materialized
end

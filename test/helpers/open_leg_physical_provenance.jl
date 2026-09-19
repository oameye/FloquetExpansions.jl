if !isdefined(@__MODULE__, :OpenLegEvaluationPlan)
  include(joinpath(@__DIR__, "open_leg_evaluation_plan.jl"))
end
if !isdefined(@__MODULE__, :StinespringHistory)
  include(joinpath(@__DIR__, "stinespring_normalization_reference.jl"))
end

struct PhysicalOpenLegVertex{C,H,S,T}
  order::Int
  channel::C
  system_harmonic::H
  output_sideband::S
  mismatch::H
  grade::Int
  value::T
end

function physical_jump_vertex(
  order::Int, channel::C, system_harmonic::H, output_sideband::H, value::T
) where {C,H,T}
  return PhysicalOpenLegVertex(
    order,
    channel,
    system_harmonic,
    output_sideband,
    system_harmonic - output_sideband,
    1,
    value,
  )
end

function physical_drift_vertex(
  order::Int, channel::C, system_harmonic::H, output_sideband::H, value::T
) where {C,H,T}
  return PhysicalOpenLegVertex(
    order, channel, system_harmonic, output_sideband, system_harmonic, 0, value
  )
end

struct PhysicalOpenLegCompiled{H,T,V}
  input::OpenLegCompiledInput{H,T}
  provenance::Vector{V}
end

function compile_physical_openleg_input(
  vertices_by_order::Vector{Vector{V}}, zero_component::T
) where {C,H,T,V<:PhysicalOpenLegVertex{C,H,H,T}}
  paths = Tuple[()]
  path_ids = Dict{Tuple,Int}(() => 1)
  values = T[]
  provenance = V[]
  compiled_by_order = Vector{Vector{OpenLegPlanVertex{H}}}()

  for (order, vertices) in enumerate(vertices_by_order)
    compiled = OpenLegPlanVertex{H}[]
    for vertex in vertices
      vertex.order == order || throw(ArgumentError("vertex order does not match its tier"))
      vertex.grade in (0, 1) ||
        throw(ArgumentError("bare physical reference vertices must have grade zero or one"))
      legs = if iszero(vertex.grade)
        ()
      else
        (OpenLegOutputLeg(vertex.channel, vertex.output_sideband),)
      end
      path = openleg_intern_path!(paths, path_ids, legs)
      push!(values, vertex.value)
      push!(provenance, vertex)
      push!(
        compiled,
        OpenLegPlanVertex(order, vertex.mismatch, path, vertex.grade, length(values)),
      )
    end
    push!(compiled_by_order, compiled)
  end

  input = OpenLegCompiledInput(compiled_by_order, values, paths, zero_component)
  return PhysicalOpenLegCompiled(input, provenance)
end

function compile_physical_openleg_evaluation_plan(
  vertices_by_order::Vector{Vector{V}}, order::Int, zero_component::T
) where {C,T,V<:PhysicalOpenLegVertex{C,Int,Int,T}}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  physical = compile_physical_openleg_input(vertices_by_order, zero_component)
  input = physical.input
  zero_harmonic = 0
  paths = input.paths
  path_ids = Dict{Tuple,Int}(path => index for (index, path) in enumerate(paths))
  identity_state = OpenLegPlanState(zero_harmonic, 1, 0)

  wave_support = [OpenLegPlanState{Int}[] for _ in 1:order]
  effective_support = [OpenLegPlanState{Int}[] for _ in 1:order]
  residuals = [OpenLegResidualNode{Int}[] for _ in 1:order]

  residual_node_count = 0
  generator_product_count = 0
  fold_product_count = 0

  for n in 1:order
    generator_edges = Dict{OpenLegPlanState{Int},Vector{OpenLegGeneratorEdge{Int}}}()
    fold_edges = Dict{OpenLegPlanState{Int},Vector{OpenLegFoldEdge{Int}}}()
    residual_support = OpenLegPlanState{Int}[]

    for vertex_order in 1:min(n, length(input.vertices_by_order))
      previous_order = n - vertex_order
      previous_states = if iszero(previous_order)
        OpenLegPlanState{Int}[identity_state]
      else
        wave_support[previous_order]
      end
      for vertex in input.vertices_by_order[vertex_order], previous in previous_states
        output = openleg_compose_state!(
          vertex.harmonic, vertex.path, vertex.grade, previous, paths, path_ids
        )
        openleg_push_unique_state!(residual_support, output)
        edges = get!(generator_edges, output) do
          return OpenLegGeneratorEdge{Int}[]
        end
        push!(edges, OpenLegGeneratorEdge(vertex.component, previous_order, previous))
        generator_product_count += 1
      end
    end

    for wave_order in 1:(n - 1)
      effective_order = n - wave_order
      for wave_state in wave_support[wave_order],
        effective_state in effective_support[effective_order]

        output = openleg_compose_states!(wave_state, effective_state, paths, path_ids)
        openleg_push_unique_state!(residual_support, output)
        edges = get!(fold_edges, output) do
          return OpenLegFoldEdge{Int}[]
        end
        push!(
          edges, OpenLegFoldEdge(wave_order, wave_state, effective_order, effective_state)
        )
        fold_product_count += 1
      end
    end

    for state in residual_support
      push!(
        residuals[n],
        OpenLegResidualNode(
          n,
          state,
          get(generator_edges, state, OpenLegGeneratorEdge{Int}[]),
          get(fold_edges, state, OpenLegFoldEdge{Int}[]),
        ),
      )
      residual_node_count += 1
    end

    effective_support[n] = OpenLegPlanState{Int}[
      state for state in residual_support if openleg_is_model_state(state)
    ]
    wave_support[n] = OpenLegPlanState{Int}[
      state for state in residual_support if !openleg_is_model_state(state)
    ]
  end

  compiled_input = OpenLegCompiledInput(
    input.vertices_by_order, input.values, paths, input.zero_component
  )
  counts = OpenLegEvaluationPlanCounts(
    residual_node_count,
    sum(length, wave_support),
    sum(length, effective_support),
    generator_product_count,
    fold_product_count,
    length(paths),
  )
  plan = OpenLegEvaluationPlan(
    compiled_input, order, zero_harmonic, wave_support, effective_support, residuals, counts
  )
  return (; plan, provenance=physical.provenance)
end

function physical_openleg_output_key(path::Tuple)
  isempty(path) && return ()
  ordered = sort(collect(path); by=leg -> (string(leg.channel), leg.sideband))
  return Tuple(ordered)
end

struct PhysicalOpenLegGrowthCounts
  raw_contributions::Int
  history_paths::Int
  physical_output_classes::Int
  recurrence_states::Int
  dag_nodes::Int
  generator_products::Int
  fold_products::Int
  physical_output_coalescences::Int
  left_right_paste_pairings::Int
end

function physical_openleg_growth_counts(plan::OpenLegEvaluationPlan)
  class_multiplicity = Dict{Tuple,Int}()
  history_paths = 0
  for path in plan.input.paths
    isempty(path) && continue
    history_paths += 1
    key = physical_openleg_output_key(path)
    class_multiplicity[key] = get(class_multiplicity, key, 0) + 1
  end

  output_classes = length(class_multiplicity)
  coalescences = sum(count - 1 for count in values(class_multiplicity))
  paste_pairings = sum(count^2 for count in values(class_multiplicity))
  recurrence_states = plan.counts.wave_nodes + plan.counts.effective_nodes
  raw_contributions = plan.counts.generator_products + plan.counts.fold_products
  return PhysicalOpenLegGrowthCounts(
    raw_contributions,
    history_paths,
    output_classes,
    recurrence_states,
    plan.counts.residual_nodes,
    plan.counts.generator_products,
    plan.counts.fold_products,
    coalescences,
    paste_pairings,
  )
end

function physical_openleg_collapse_vertices(vertices_by_order, zero_component)
  amplitudes = OpenLegPeriodic[]
  for vertices in vertices_by_order
    components = Dict{Tuple{Int,Int},typeof(zero_component)}()
    for vertex in vertices
      key = (vertex.mismatch, vertex.grade)
      components[key] = get(components, key, zero_component) + vertex.value
    end
    push!(amplitudes, openleg_periodic(components, zero_component))
  end
  return amplitudes
end

function halfchain_rr_pasted_component(
  amplitudes::AbstractDict{Int,T},
  density_harmonic::Int;
  product,
  paste,
  minus_imaginary,
  zero_component,
) where {T}
  density_harmonic > 0 || throw(ArgumentError("density harmonic must be positive"))
  result = zero_component
  for (p, Lp) in amplitudes, (q, Lq) in amplitudes
    p - q == density_harmonic || continue
    for (r, Lr) in amplitudes, (s, Ls) in amplitudes
      r - s == -density_harmonic || continue
      forward = paste(product(Lp, Lr), product(Lq, Ls))
      reverse = paste(product(Lr, Lp), product(Ls, Lq))
      result += minus_imaginary * (1 // density_harmonic) * (forward - reverse)
    end
  end
  return result
end

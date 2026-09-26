struct StaticExpDependency{H}
  left_order::Int
  left_harmonic::H
  right_order::Int
  right_harmonic::H
  right_node::Int
end

struct StaticExpPlan{H}
  zero_harmonic::H
  nodes::Vector{Vector{StaticExpDependency{H}}}
  targets::Vector{Vector{Int}}
end

struct VanVleckNormalizationPlan{H,C}
  connected_log::ConnectedLogPlan{H,C}
  static_factor::StaticExpPlan{H}
end

function static_power_support(plan::ConnectedLogPlan{H}) where {H}
  highest_power = length(plan.outputs)
  supports = [[Set{H}() for _ in 1:highest_power] for _ in 1:highest_power]
  for n in 1:highest_power
    union!(supports[1][n], keys(plan.outputs[n]))
  end
  for power in 2:highest_power, n in power:highest_power, k in 1:(n - power + 1)
    for left in supports[1][k], right in supports[power - 1][n - k]
      push!(supports[power][n], left + right)
    end
  end
  return supports
end

function static_exp_node!(
  state::Tuple{Int,Int,H},
  supports::Vector{Vector{Set{H}}},
  nodes::Vector{Vector{StaticExpDependency{H}}},
  node_ids::Dict{Tuple{Int,Int,H},Int},
) where {H}
  haskey(node_ids, state) && return node_ids[state]
  power, n, harmonic = state
  dependencies = StaticExpDependency{H}[]
  for k in 1:(n - power + 1), left_harmonic in supports[1][k]
    right_harmonic = harmonic - left_harmonic
    right_harmonic in supports[power - 1][n - k] || continue
    right_node = if power == 2
      0
    else
      static_exp_node!((power - 1, n - k, right_harmonic), supports, nodes, node_ids)
    end
    push!(
      dependencies, StaticExpDependency(k, left_harmonic, n - k, right_harmonic, right_node)
    )
  end
  push!(nodes, dependencies)
  node_ids[state] = length(nodes)
  return length(nodes)
end

function compile_static_exp_plan(plan::ConnectedLogPlan{H}, zero_harmonic::H) where {H}
  highest_power = length(plan.outputs)
  supports = static_power_support(plan)
  nodes = Vector{StaticExpDependency{H}}[]
  node_ids = Dict{Tuple{Int,Int,H},Int}()
  targets = [zeros(Int, highest_power) for _ in 1:highest_power]
  for n in 2:highest_power, power in 2:n
    zero_harmonic in supports[power][n] || continue
    targets[n][power] = static_exp_node!(
      (power, n, zero_harmonic), supports, nodes, node_ids
    )
  end
  return StaticExpPlan(zero_harmonic, nodes, targets)
end

function compile_van_vleck_normalization_plan(plan::WaveOperatorPlan)
  connected_log = compile_connected_log_plan(plan)
  static_factor = compile_static_exp_plan(connected_log, plan.zero_harmonic)
  return VanVleckNormalizationPlan(connected_log, static_factor)
end

function evaluate_static_factor(
  plan::StaticExpPlan{H},
  connected_log::Vector{Dict{H,T}},
  identity_component::T,
  zero_component::T,
  product,
) where {H,T}
  values = Vector{T}(undef, length(plan.nodes))
  for (node, dependencies) in enumerate(plan.nodes)
    value = zero_component
    for dependency in dependencies
      left = get(
        connected_log[dependency.left_order], dependency.left_harmonic, zero_component
      )
      right = if dependency.right_node == 0
        get(connected_log[dependency.right_order], dependency.right_harmonic, zero_component)
      else
        values[dependency.right_node]
      end
      value += product(left, right)
    end
    values[node] = simplify_component(value)::T
  end

  static_factor = T[identity_component]
  for n in eachindex(connected_log)
    value = get(connected_log[n], plan.zero_harmonic, zero_component)
    for power in 2:n
      node = plan.targets[n][power]
      node == 0 && continue
      value += (1 // factorial(power)) * values[node]
    end
    push!(static_factor, simplify_component(value)::T)
  end
  return static_factor
end

function static_series_inverse(series::Vector{T}, highest_power::Int, product) where {T}
  inverse = T[first(series)]
  for n in 1:highest_power
    coefficient = zero(first(series))
    for k in 1:n
      coefficient += product(series[k + 1], inverse[n - k + 1])
    end
    push!(inverse, simplify_component(-coefficient)::T)
  end
  return inverse
end

function static_series_product(
  left::Vector{T}, right::Vector{T}, highest_power::Int, product
) where {T}
  result = Vector{T}(undef, highest_power + 1)
  for n in 0:highest_power
    coefficient = zero(first(left))
    for k in 0:n
      coefficient += product(left[k + 1], right[n - k + 1])
    end
    result[n + 1] = simplify_component(coefficient)::T
  end
  return result
end

function normalize_to_van_vleck(
  plan::VanVleckNormalizationPlan{H},
  wave_operator::WaveOperator{H,T},
  components::AbstractDict{H,T},
  zero_component::T,
  conventions::BlochConventions,
) where {H,T}
  (; product, weight_phase) = conventions
  bloch_effective = wave_operator.bloch_effective_generator
  highest_power = length(bloch_effective) - 1

  connected_log = evaluate_connected_log(
    plan.connected_log, components, zero_component, product
  )
  apply_weight_phase!(connected_log, weight_phase)

  static_factor = evaluate_static_factor(
    plan.static_factor, connected_log, one(first(bloch_effective)), zero_component, product
  )
  inverse_static_factor = static_series_inverse(static_factor, highest_power, product)
  bloch_effective_times_static_factor = static_series_product(
    bloch_effective, static_factor, highest_power, product
  )
  effective = static_series_product(
    inverse_static_factor, bloch_effective_times_static_factor, highest_power, product
  )
  return (; connected_log, effective)
end

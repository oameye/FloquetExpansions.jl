using FloquetExpansions

isdefined(@__MODULE__, :LyndonLogEvaluationPlan) ||
  include(joinpath(@__DIR__, "lyndon_log_evaluator.jl"))

struct StaticExpDependency{H}
  left_order::Int
  left_harmonic::H
  right_order::Int
  right_harmonic::H
  right_node::Int
end

struct StaticExpNode{H}
  power::Int
  order::Int
  harmonic::H
  dependencies::Vector{StaticExpDependency{H}}
end

struct StaticSectorExpPlan{H}
  zero_harmonic::H
  nodes::Vector{StaticExpNode{H}}
  targets::Vector{Vector{Int}}
  product_count::Int
end

function static_exp_power_support(log_plan::LyndonLogEvaluationPlan{H}) where {H}
  order = length(log_plan.outputs)
  supports = [[Set{H}() for _ in 1:order] for _ in 1:order]

  for n in 1:order
    union!(supports[1][n], keys(log_plan.outputs[n]))
  end

  for power in 2:order, n in power:order
    for k in 1:(n - power + 1)
      for left in supports[1][k], right in supports[power - 1][n - k]
        push!(supports[power][n], left + right)
      end
    end
  end

  return supports
end

function compile_static_sector_exp_plan(
  log_plan::LyndonLogEvaluationPlan{H}, zero_harmonic::H
) where {H}
  order = length(log_plan.outputs)
  supports = static_exp_power_support(log_plan)
  nodes = StaticExpNode{H}[]
  node_ids = Dict{Tuple{Int,Int,H},Int}()
  product_count = Ref(0)

  function ensure_node(power::Int, n::Int, harmonic::H)
    state = (power, n, harmonic)
    haskey(node_ids, state) && return node_ids[state]
    dependencies = StaticExpDependency{H}[]

    for k in 1:(n - power + 1)
      right_order = n - k
      for left_harmonic in supports[1][k]
        for right_harmonic in supports[power - 1][right_order]
          left_harmonic + right_harmonic == harmonic || continue
          right_node = power == 2 ? 0 : ensure_node(power - 1, right_order, right_harmonic)
          push!(
            dependencies,
            StaticExpDependency(k, left_harmonic, right_order, right_harmonic, right_node),
          )
          product_count[] += 1
        end
      end
    end

    push!(nodes, StaticExpNode(power, n, harmonic, dependencies))
    node = length(nodes)
    node_ids[state] = node
    return node
  end

  targets = [zeros(Int, order) for _ in 1:order]
  for n in 2:order, power in 2:n
    zero_harmonic in supports[power][n] || continue
    targets[n][power] = ensure_node(power, n, zero_harmonic)
  end

  return StaticSectorExpPlan(zero_harmonic, nodes, targets, product_count[])
end

function evaluate_static_sector_exp_plan(
  plan::StaticSectorExpPlan{H},
  log_embedding::Vector{Dict{H,T}},
  identity_component::T,
  zero_component::T;
  product,
  simplifier=identity,
) where {H,T}
  values = Vector{T}(undef, length(plan.nodes))

  for (node_id, node) in enumerate(plan.nodes)
    value = zero_component
    for dependency in node.dependencies
      left = get(
        log_embedding[dependency.left_order], dependency.left_harmonic, zero_component
      )
      right = if dependency.right_node == 0
        get(log_embedding[dependency.right_order], dependency.right_harmonic, zero_component)
      else
        values[dependency.right_node]
      end
      value += product(left, right)
    end
    values[node_id] = simplifier(value)::T
  end

  static_factor = T[identity_component]
  for n in eachindex(log_embedding)
    value = get(log_embedding[n], plan.zero_harmonic, zero_component)
    for power in 2:n
      node = plan.targets[n][power]
      node == 0 && continue
      value += (1 // factorial(power)) * values[node]
    end
    push!(static_factor, simplifier(value)::T)
  end

  return static_factor
end

using FloquetExpansions

isdefined(@__MODULE__, :LyndonLogEvaluationPlan) ||
  include(joinpath(@__DIR__, "lyndon_log_evaluator.jl"))

const FE_CVVR = FloquetExpansions

struct ConnectedVanVleckResult{H,T}
  static_factor::Vector{T}
  inverse_static_factor::Vector{T}
  log_embedding::Vector{Dict{H,T}}
  effective::Vector{T}
  counts::FE_CVVR.BlochVanVleckCounts
end

function connected_static_factor(
  log_embedding::Vector{Dict{H,T}},
  zero_harmonic::H,
  identity_component::T,
  zero_component::T;
  product,
  simplifier=identity,
) where {H,T}
  order = length(log_embedding)
  powers = [[Dict{H,T}() for _ in 1:max(order, 1)] for _ in 1:max(order, 1)]
  static_factor = T[identity_component]
  counts = FE_CVVR.BlochVanVleckCounts()

  for n in 1:order
    powers[1][n] = log_embedding[n]
    static_n = get(log_embedding[n], zero_harmonic, zero_component)

    for power in 2:n
      power_coefficient = Dict{H,T}()
      for k in 1:(n - power + 1)
        counts.log_products += 1
        contribution = FE_CVVR.bloch_vv_periodic_product(
          log_embedding[k], powers[power - 1][n - k], product, simplifier, counts
        )
        power_coefficient = FE_CVVR.bloch_vv_add(
          power_coefficient, contribution, simplifier
        )
      end
      powers[power][n] = power_coefficient
      weight = 1 // factorial(power)
      static_n = simplifier(
        static_n + weight * get(power_coefficient, zero_harmonic, zero_component)
      )::T
    end

    push!(static_factor, simplifier(static_n)::T)
  end

  return static_factor, counts
end

function connected_van_vleck_reconstruction(
  projection_plan::FE_CVVR.BlochProjectionPlan{H},
  bloch::FE_CVVR.BlochProjectionResult{H,T},
  components::AbstractDict{H,T},
  log_plan::LyndonLogEvaluationPlan{H};
  product,
  zero_component::T,
  simplifier=identity,
) where {H,T}
  order = length(bloch.effective) - 1
  projection_plan.order == order + 1 ||
    throw(ArgumentError("Bloch plan/result truncations are inconsistent"))
  length(log_plan.outputs) == order ||
    throw(ArgumentError("connected-log plan/result truncations are inconsistent"))

  log_embedding = evaluate_lyndon_log_plan(
    log_plan, components; product, zero_component, simplifier
  )
  identity_component = one(first(bloch.effective))
  static_factor, counts = connected_static_factor(
    log_embedding,
    projection_plan.zero_harmonic,
    identity_component,
    zero_component;
    product,
    simplifier,
  )
  inverse_static_factor = FE_CVVR.bloch_vv_static_series_inverse(
    static_factor, order, product, counts; simplifier
  )
  right_transformed = FE_CVVR.bloch_vv_static_series_product(
    bloch.effective, static_factor, order, product, counts; simplifier
  )
  effective = FE_CVVR.bloch_vv_static_series_product(
    inverse_static_factor, right_transformed, order, product, counts; simplifier
  )

  return ConnectedVanVleckResult(
    static_factor, inverse_static_factor, log_embedding, effective, counts
  )
end

function connected_reconstruction_backend_products(
  result::ConnectedVanVleckResult, log_plan::LyndonLogEvaluationPlan
)
  counts = result.counts
  return lyndon_backend_products(log_plan) +
         counts.harmonic_products +
         counts.inverse_products +
         counts.similarity_products
end

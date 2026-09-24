struct CKHoriDepritResult{K}
  generator::Vector{K}
  derivative::Vector{K}
  effective::Vector{K}
  products::Int
end

function ck_endpoint_solve_homological(state::CKEndpointKernel{H,T}) where {H,T}
  result = Dict{CKEndpointKey{H},T}()
  for (key, value) in state.terms
    isempty(key.vertices) && continue
    whole = CKEndpointConstraint(0, length(key.vertices))
    whole in key.model_intervals && continue

    resolvent_intervals = copy(key.resolvent_intervals)
    push!(resolvent_intervals, whole)
    solved = CKEndpointKey(
      copy(key.vertices), copy(key.model_intervals), resolvent_intervals
    )
    ck_endpoint_accumulate!(result, solved, value, state.zero_component)
  end
  return CKEndpointKernel(result, state.zero_component)
end

function ck_endpoint_ad_coefficient(
  generator::AbstractVector{K}, previous::AbstractVector{K}, n::Int, zero_state::K
) where {K<:CKEndpointKernel}
  coefficient = zero_state
  products = 0
  for j in 1:(n - 1)
    left = generator[j]
    right = previous[n - j]
    (isempty(left.terms) || isempty(right.terms)) && continue
    coefficient += ck_endpoint_product(left, right)
    coefficient -= ck_endpoint_product(right, left)
    products += 2
  end
  return coefficient, products
end

function ck_endpoint_hori_deprit_source!(
  amplitude_tower::Vector{Vector{K}},
  derivative_tower::Vector{Vector{K}},
  generator::Vector{K},
  n::Int,
  zero_state::K,
) where {K<:CKEndpointKernel}
  source = amplitude_tower[1][n]
  products = 0

  for power in 1:(n - 1)
    coefficient, count = ck_endpoint_ad_coefficient(
      generator, amplitude_tower[power], n, zero_state
    )
    amplitude_tower[power + 1][n] = coefficient
    source += ck_endpoint_scale(coefficient, (-1)^power // factorial(power))
    products += count
  end

  for power in 1:(n - 1)
    coefficient, count = ck_endpoint_ad_coefficient(
      generator, derivative_tower[power], n, zero_state
    )
    derivative_tower[power + 1][n] = coefficient
    weight = (-1)^(power + 1) // factorial(power + 1)
    source += ck_endpoint_scale(coefficient, weight)
    products += count
  end

  return source, products
end

function evaluate_ck_endpoint_hori_deprit(
  amplitude_orders::AbstractVector{K}, order::Int, zero_state::K
) where {K<:CKEndpointKernel}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  isempty(amplitude_orders) && throw(ArgumentError("amplitude series must not be empty"))

  amplitude_tower = [K[zero_state for _ in 1:order] for _ in 1:order]
  derivative_tower = [K[zero_state for _ in 1:order] for _ in 1:order]
  for n in 1:min(order, length(amplitude_orders))
    amplitude_tower[1][n] = amplitude_orders[n]
  end

  generator = K[zero_state for _ in 1:order]
  derivative = derivative_tower[1]
  effective = K[zero_state for _ in 1:order]
  products = 0

  for n in 1:order
    source, count = ck_endpoint_hori_deprit_source!(
      amplitude_tower, derivative_tower, generator, n, zero_state
    )
    products += count
    effective[n] = ck_endpoint_project_model(source)
    derivative[n] = source - effective[n]
    generator[n] = ck_endpoint_solve_homological(source)
  end

  return CKHoriDepritResult(generator, derivative, effective, products)
end

function evaluate_ck_hori_deprit(
  amplitude_orders::AbstractVector{K}, order::Int, zero_state::K
) where {K<:CKPhysicalKernel}
  isempty(amplitude_orders) && throw(ArgumentError("amplitude series must not be empty"))
  endpoint_orders = [ck_endpoint_kernel(amplitude) for amplitude in amplitude_orders]
  zero_endpoint = ck_endpoint_kernel(zero_state)
  return evaluate_ck_endpoint_hori_deprit(endpoint_orders, order, zero_endpoint)
end

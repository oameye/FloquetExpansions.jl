using LinearAlgebra: kron

struct OpenLegPeriodic{T}
  components::Dict{Tuple{Int,Int},T}
  zero_component::T
end

openleg_component_iszero(value) = iszero(value)
openleg_component_iszero(value::AbstractMatrix) = all(iszero, value)

function openleg_periodic(components::Dict{Tuple{Int,Int},T}, zero_component::T) where {T}
  cleaned = Dict{Tuple{Int,Int},T}()
  for (key, value) in components
    openleg_component_iszero(value) || (cleaned[key] = value)
  end
  return OpenLegPeriodic(cleaned, zero_component)
end

Base.zero(value::OpenLegPeriodic{T}) where {T} =
  OpenLegPeriodic(Dict{Tuple{Int,Int},T}(), zero(value.zero_component))

function Base.getindex(value::OpenLegPeriodic, harmonic::Int, grade::Int)
  return get(value.components, (harmonic, grade), value.zero_component)
end

function Base.:+(left::OpenLegPeriodic{T}, right::OpenLegPeriodic{T}) where {T}
  out = copy(left.components)
  for (key, value) in right.components
    out[key] = get(out, key, left.zero_component) + value
  end
  return openleg_periodic(out, left.zero_component)
end

function Base.:-(left::OpenLegPeriodic{T}, right::OpenLegPeriodic{T}) where {T}
  out = copy(left.components)
  for (key, value) in right.components
    out[key] = get(out, key, left.zero_component) - value
  end
  return openleg_periodic(out, left.zero_component)
end

function Base.:-(value::OpenLegPeriodic{T}) where {T}
  return openleg_periodic(
    Dict(key => -component for (key, component) in value.components), value.zero_component
  )
end

function Base.:*(scalar::Number, value::OpenLegPeriodic{T}) where {T}
  return openleg_periodic(
    Dict(key => scalar * component for (key, component) in value.components),
    value.zero_component,
  )
end

Base.:*(value::OpenLegPeriodic, scalar::Number) = scalar * value

function openleg_identity(template::OpenLegPeriodic{T}, identity_component::T) where {T}
  return openleg_periodic(Dict((0, 0) => identity_component), template.zero_component)
end

function openleg_product(left::OpenLegPeriodic{T}, right::OpenLegPeriodic{T}, product) where {T}
  out = Dict{Tuple{Int,Int},T}()
  for ((left_harmonic, left_grade), left_component) in left.components
    for ((right_harmonic, right_grade), right_component) in right.components
      key = (left_harmonic + right_harmonic, left_grade + right_grade)
      term = product(left_component, right_component)
      out[key] = get(out, key, left.zero_component) + term
    end
  end
  return openleg_periodic(out, left.zero_component)
end

function openleg_commutator(left::OpenLegPeriodic, right::OpenLegPeriodic, product)
  return openleg_product(left, right, product) - openleg_product(right, left, product)
end

function openleg_project(value::OpenLegPeriodic{T}) where {T}
  return openleg_periodic(
    Dict(key => component for (key, component) in value.components if iszero(first(key))),
    value.zero_component,
  )
end

function openleg_q_inverse(value::OpenLegPeriodic{T}) where {T}
  out = Dict{Tuple{Int,Int},T}()
  for ((harmonic, grade), component) in value.components
    iszero(harmonic) && continue
    out[(harmonic, grade)] = (im * (1 // harmonic)) * component
  end
  return openleg_periodic(out, value.zero_component)
end

function openleg_derivative(value::OpenLegPeriodic{T}) where {T}
  out = Dict{Tuple{Int,Int},T}()
  for ((harmonic, grade), component) in value.components
    iszero(harmonic) && continue
    out[(harmonic, grade)] = (-im * harmonic) * component
  end
  return openleg_periodic(out, value.zero_component)
end

function openleg_equal(left::OpenLegPeriodic, right::OpenLegPeriodic)
  all_keys = union(keys(left.components), keys(right.components))
  return all(left[key...] == right[key...] for key in all_keys)
end

function openleg_grades(value::OpenLegPeriodic)
  return sort!(unique(last(key) for key in keys(value.components)))
end

struct OpenLegBlochReference{P}
  wave::Vector{P}
  effective::Vector{P}
end

function openleg_bloch_reference(
  amplitude_orders::Vector{P}, order::Int; product, identity_component
) where {P<:OpenLegPeriodic}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  isempty(amplitude_orders) && throw(ArgumentError("at least one amplitude order is required"))

  template = first(amplitude_orders)
  identity = openleg_identity(template, identity_component)
  wave = P[]
  effective = P[]

  for n in 1:order
    residual = zero(template)
    for r in 1:min(n, length(amplitude_orders))
      previous = n == r ? identity : wave[n - r]
      residual = residual + openleg_product(amplitude_orders[r], previous, product)
    end
    for j in 1:(n - 1)
      residual = residual - openleg_product(wave[j], effective[n - j], product)
    end

    push!(effective, openleg_project(residual))
    push!(wave, openleg_q_inverse(residual))
  end

  return OpenLegBlochReference(wave, effective)
end

struct OpenLegHoriDepritOrder2{P}
  generator1::P
  generator2::P
  effective1::P
  effective2::P
end

function openleg_hori_deprit_order2(A1::P, A2::P; product) where {P<:OpenLegPeriodic}
  B1 = openleg_project(A1)
  G1 = openleg_q_inverse(A1)
  derivative_G1 = openleg_derivative(G1)
  F2 =
    A2 -
    openleg_commutator(G1, A1, product) +
    (1 // 2) * openleg_commutator(G1, derivative_G1, product)
  B2 = openleg_project(F2)
  G2 = openleg_q_inverse(F2)
  return OpenLegHoriDepritOrder2(G1, G2, B1, B2)
end

function static_unitary_jump_tp_coefficient(order::Int)
  order >= 0 || throw(ArgumentError("order must be nonnegative"))
  result = 0 // 1
  for q in 0:order
    for left_order in 0:(order - q)
      right_order = order - q - left_order
      result +=
        (1 // factorial(q)) *
        ((-1 // 2)^left_order / factorial(left_order)) *
        ((-1 // 2)^right_order / factorial(right_order))
    end
  end
  return result
end

function openleg_triangle_contains_channel_order(generator_order::Int, channel_order::Int)
  generator_order >= 0 || throw(ArgumentError("generator order must be nonnegative"))
  channel_order >= 0 || throw(ArgumentError("channel order must be nonnegative"))
  channel_order <= generator_order + 1 || return false

  for q in 0:channel_order
    for left_order in 0:(channel_order - q)
      right_order = channel_order - q - left_order
      q + left_order <= generator_order + 1 || return false
      q + right_order <= generator_order + 1 || return false
    end
  end
  return true
end

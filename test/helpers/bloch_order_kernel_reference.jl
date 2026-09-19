struct BlochOrderResult{T}
  wave::Vector{T}
  effective::Vector{T}
  products::Int
end

function bloch_order_recurrence(
  generators_by_order::Vector{T},
  order::Int,
  identity_component::T,
  zero_component::T;
  product,
  project_model,
  solve_complement,
  simplifier=identity,
) where {T}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  wave = T[]
  effective = T[]
  products = 0

  for n in 1:order
    residual = zero_component

    for generator_order in 1:min(n, length(generators_by_order))
      previous_order = n - generator_order
      previous = iszero(previous_order) ? identity_component : wave[previous_order]
      residual += product(generators_by_order[generator_order], previous)
      products += 1
    end

    for wave_order in 1:(n - 1)
      effective_order = n - wave_order
      residual -= product(wave[wave_order], effective[effective_order])
      products += 1
    end

    residual = simplifier(residual)::T
    push!(effective, simplifier(project_model(residual))::T)
    if n < order
      push!(wave, simplifier(solve_complement(residual))::T)
    end
  end

  return BlochOrderResult(wave, effective, products)
end

struct KernelFourierSeries{T}
  components::Dict{Int,T}
  zero_component::T
end

function kernel_fourier_series(components::AbstractDict{Int,T}, zero_component::T) where {T}
  return KernelFourierSeries(
    Dict(harmonic => value for (harmonic, value) in components if value != zero_component),
    zero_component,
  )
end

function Base.:+(left::KernelFourierSeries{T}, right::KernelFourierSeries{T}) where {T}
  result = copy(left.components)
  for (harmonic, value) in right.components
    updated = get(result, harmonic, left.zero_component) + value
    if updated == left.zero_component
      haskey(result, harmonic) && delete!(result, harmonic)
    else
      result[harmonic] = updated
    end
  end
  return KernelFourierSeries(result, left.zero_component)
end

function Base.:-(left::KernelFourierSeries{T}, right::KernelFourierSeries{T}) where {T}
  result = copy(left.components)
  for (harmonic, value) in right.components
    updated = get(result, harmonic, left.zero_component) - value
    if updated == left.zero_component
      haskey(result, harmonic) && delete!(result, harmonic)
    else
      result[harmonic] = updated
    end
  end
  return KernelFourierSeries(result, left.zero_component)
end

function kernel_fourier_product(
  left::KernelFourierSeries{T}, right::KernelFourierSeries{T}
) where {T}
  result = Dict{Int,T}()
  for (left_harmonic, left_value) in left.components,
    (right_harmonic, right_value) in right.components

    harmonic = left_harmonic + right_harmonic
    updated = get(result, harmonic, left.zero_component) + left_value * right_value
    if updated == left.zero_component
      haskey(result, harmonic) && delete!(result, harmonic)
    else
      result[harmonic] = updated
    end
  end
  return KernelFourierSeries(result, left.zero_component)
end

function kernel_fourier_project_model(series::KernelFourierSeries{T}) where {T}
  value = get(series.components, 0, series.zero_component)
  return value == series.zero_component ?
         KernelFourierSeries(Dict{Int,T}(), series.zero_component) :
         KernelFourierSeries(Dict(0 => value), series.zero_component)
end

function kernel_fourier_solve_complement(series::KernelFourierSeries{T}) where {T}
  result = Dict{Int,T}()
  for (harmonic, value) in series.components
    iszero(harmonic) && continue
    result[harmonic] = (im // harmonic) * value
  end
  return KernelFourierSeries(result, series.zero_component)
end

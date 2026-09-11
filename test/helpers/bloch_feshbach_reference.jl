using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra

mutable struct BlochReferenceCounts
  full_products::Int
  static_products::Int
  component_products::Int
end

BlochReferenceCounts() = BlochReferenceCounts(0, 0, 0)

struct BlochReferenceResult{P,E}
  wave::Vector{P}
  effective::Vector{E}
  counts::BlochReferenceCounts
end

function reference_periodic_product(
  L::PeriodicGenerator{T}, R::PeriodicGenerator{T}, product, counts::BlochReferenceCounts
) where {T}
  isequal(L.wd, R.wd) ||
    throw(ArgumentError("periodic generators use different frequencies"))

  out = Dict{Int,T}()
  for left_harmonic in keys(L), right_harmonic in keys(R)
    counts.component_products += 1
    harmonic = left_harmonic + right_harmonic
    term = product(L[left_harmonic], R[right_harmonic])
    out[harmonic] = haskey(out, harmonic) ? out[harmonic] + term : term
  end
  return PeriodicGenerator(out, L.wd, L.zero_component)
end

function reference_right_static_product(
  X::PeriodicGenerator{T}, B::T, product, counts::BlochReferenceCounts
) where {T}
  out = Dict{Int,T}()
  for harmonic in keys(X)
    counts.component_products += 1
    out[harmonic] = product(X[harmonic], B)
  end
  return PeriodicGenerator(out, X.wd, X.zero_component)
end

function reference_q_inverse(G::PeriodicGenerator{T}, inverse_weight) where {T}
  out = Dict{Int,T}()
  for harmonic in keys(G)
    iszero(harmonic) && continue
    out[harmonic] = inverse_weight(harmonic) * G[harmonic]
  end
  return PeriodicGenerator(out, G.wd, G.zero_component)
end

function bloch_reference(
  G::PeriodicGenerator{T}, order::Int; product, inverse_weight
) where {T}
  order >= 1 || throw(ArgumentError("order must be >= 1"))

  wave = PeriodicGenerator{T}[]
  effective = T[]
  counts = BlochReferenceCounts()

  residual = G
  push!(effective, SQA.simplify(time_average(residual)))
  if order > 1
    push!(wave, SQA.simplify(reference_q_inverse(residual, inverse_weight)))
  end

  for n in 1:(order - 1)
    counts.full_products += 1
    residual = reference_periodic_product(G, wave[n], product, counts)

    for j in 1:n
      counts.static_products += 1
      correction = reference_right_static_product(
        wave[j], effective[n - j + 1], product, counts
      )
      residual = residual - correction
    end

    residual = SQA.simplify(residual)
    push!(effective, SQA.simplify(time_average(residual)))

    if n < order - 1
      push!(wave, SQA.simplify(reference_q_inverse(residual, inverse_weight)))
    end
  end

  return BlochReferenceResult(wave, effective, counts)
end

function reference_static_generator(component::T, template::PeriodicGenerator{T}) where {T}
  return PeriodicGenerator(Dict(0 => component), template.wd, template.zero_component)
end

function reference_wave_equation_rhs(
  result::BlochReferenceResult{P,T}, template::P, n::Int, product
) where {T,P<:PeriodicGenerator{T}}
  rhs = reference_static_generator(result.effective[n + 1], template)
  counts = BlochReferenceCounts()
  for j in 1:n
    rhs =
      rhs + reference_right_static_product(
        result.wave[j], result.effective[n - j + 1], product, counts
      )
  end
  return SQA.simplify(rhs)
end

function expected_bloch_series_products(order::Int)
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  return (order - 1) * (order + 2) ÷ 2
end

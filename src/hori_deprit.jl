triindex(n::Int, j::Int) = (n * (n + 1)) ÷ 2 + j + 1
function lie_transform_phase(generator::PeriodicGenerator)
  return -component_convention(generator).generator_phase
end

function weight_generator(generator::PeriodicGenerator, j::Int)
  return lie_transform_phase(generator)^j * (1 // factorial(j))
end

function weight_micromotion_derivative(generator::PeriodicGenerator, j::Int)
  return lie_transform_phase(generator)^j * (1 // factorial(j + 1))
end

function dressed_generator_node(
  K::Vector{P}, dressed_generator::Vector{P}, n::Int, j::Int, generator::P
) where {P}
  result = zero(generator)
  for k in 1:(n - j + 1)
    previous = dressed_generator[triindex(n - k, j - 1)]
    result = result + SQA.commutator(K[k], previous)
  end
  return result
end

function dressed_micromotion_derivative_node(
  K::Vector{P},
  Kdot::Vector{P},
  dressed_micromotion_derivative::Vector{P},
  n::Int,
  j::Int,
  generator::P,
) where {P}
  result = zero(generator)
  for k in 1:(n - j + 1)
    previous = if j == 1
      Kdot[n - k + 1]
    else
      dressed_micromotion_derivative[triindex(n - k, j - 1)]
    end
    result = result + SQA.commutator(K[k], previous)
  end
  return result
end

function assemble_resolvent(
  dressed_generator::Vector{P},
  dressed_micromotion_derivative::Vector{P},
  n::Int,
  generator::P,
) where {P<:PeriodicGenerator}
  result = zero(generator)
  for j in 0:n
    result = result + weight_generator(generator, j) * dressed_generator[triindex(n, j)]
  end
  for j in 1:n
    weight = weight_micromotion_derivative(generator, j)
    result = result - weight * dressed_micromotion_derivative[triindex(n, j)]
  end
  return result
end

function van_vleck_expansion(
  generator::P, gauge::VanVleck{HoriDeprit}, order::Int, provenance::R
) where {P<:PeriodicGenerator,R<:FloquetProvenance}
  nodes = (order * (order + 1)) ÷ 2
  dressed_generator = [zero(generator) for _ in 1:nodes]
  dressed_micromotion_derivative = [zero(generator) for _ in 1:nodes]
  generator_type = typeof(generator)
  K = generator_type[]
  Kdot = generator_type[]
  E = typeof(time_average(generator))
  effective = E[]

  for n in 0:(order - 1)
    dressed_generator[triindex(n, 0)] = n == 0 ? generator : zero(generator)

    for j in 1:n
      dressed_generator[triindex(n, j)] = dressed_generator_node(
        K, dressed_generator, n, j, generator
      )
    end

    for j in 1:n
      dressed_micromotion_derivative[triindex(n, j)] = dressed_micromotion_derivative_node(
        K, Kdot, dressed_micromotion_derivative, n, j, generator
      )
    end

    resolvent = SQA.simplify(
      assemble_resolvent(dressed_generator, dressed_micromotion_derivative, n, generator)
    )
    effective_n = SQA.simplify(time_average(resolvent))
    push!(effective, effective_n)

    if n < order - 1
      next_micromotion = SQA.simplify(antiderivative(remove_average(resolvent), gauge))
      push!(K, next_micromotion)
      push!(Kdot, derivative(next_micromotion))
    end
  end

  return FloquetExpansion(generator, K, effective, gauge, order, Uncompleted(), provenance)
end

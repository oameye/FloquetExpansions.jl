triindex(n::Int, j::Int) = (n * (n + 1)) ÷ 2 + j + 1
lie_transform_phase(::PeriodicGenerator{SQA.QAdd}) = im
lie_transform_phase(::PeriodicGenerator{Liouvillian}) = -1

function weight_generator(generator::PeriodicGenerator, j::Int)
  return lie_transform_phase(generator)^j * (1 // factorial(j))
end

function weight_kick_derivative(generator::PeriodicGenerator, j::Int)
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

function dressed_kick_derivative_node(
  K::Vector{P},
  Kdot::Vector{P},
  dressed_kick_derivative::Vector{P},
  n::Int,
  j::Int,
  generator::P,
) where {P}
  result = zero(generator)
  for k in 1:(n - j + 1)
    previous = if j == 1
      Kdot[n - k + 1]
    else
      dressed_kick_derivative[triindex(n - k, j - 1)]
    end
    result = result + SQA.commutator(K[k], previous)
  end
  return result
end

function assemble_resolvent(
  dressed_generator::Vector{P}, dressed_kick_derivative::Vector{P}, n::Int, generator::P
) where {P<:PeriodicGenerator}
  result = zero(generator)
  for j in 0:n
    result = result + weight_generator(generator, j) * dressed_generator[triindex(n, j)]
  end
  for j in 1:n
    weight = weight_kick_derivative(generator, j)
    result = result - weight * dressed_kick_derivative[triindex(n, j)]
  end
  return result
end

function floquet_expansion_impl(
  generator::P, gauge::VanVleck{HoriDeprit}, order::Int, provenance::R
) where {P<:PeriodicGenerator,R<:FloquetProvenance}
  order >= 1 || throw(ArgumentError("order must be >= 1"))

  generator isa PeriodicGenerator{SQA.QAdd} && require_hermitian_drive(generator)

  nodes = (order * (order + 1)) ÷ 2
  dressed_generator = [zero(generator) for _ in 1:nodes]
  dressed_kick_derivative = [zero(generator) for _ in 1:nodes]
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
      dressed_kick_derivative[triindex(n, j)] = dressed_kick_derivative_node(
        K, Kdot, dressed_kick_derivative, n, j, generator
      )
    end

    resolvent = SQA.simplify(
      assemble_resolvent(dressed_generator, dressed_kick_derivative, n, generator)
    )
    effective_n = SQA.simplify(time_average(resolvent))
    push!(effective, effective_n)

    if n < order - 1
      next_kick = SQA.simplify(antiderivative(remove_average(resolvent), gauge))
      push!(K, next_kick)
      push!(Kdot, derivative(next_kick))
    end
  end

  return FloquetExpansion(generator, K, effective, gauge, order, Uncompleted(), provenance)
end

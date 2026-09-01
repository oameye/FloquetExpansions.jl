using Test
using FloquetExpansions
import SecondQuantizedAlgebra as SQA
using Symbolics: Symbolics

include(joinpath(@__DIR__, "helpers", "shared.jl"))

h = FockSpace(:cavity)
a = Destroy(h, :a)
@variables w::Real g::Real ξ::Real t::Real

function drive()
  return PeriodicOperator(
    Dict(
      0 => 1 * (a' * a),
      1 => g * a,
      -1 => conj(g) * a',
      2 => ξ * (a' * a'),
      -2 => conj(ξ) * (a * a),
    ),
    w,
  )
end

comm = SQA.commutator

const SUBS = Dict(w => 1.0, g => 0.7, ξ => 0.3, t => 0.4)
function agrees(X::PeriodicOperator, Y::PeriodicOperator)
  return all(
    maxcoeff(SQA.simplify(X[l] - Y[l]), SUBS) < 1.0e-12 for l in union(keys(X), keys(Y))
  )
end

function compositions(n::Int, j::Int)
  j == 0 && return n == 0 ? [Int[]] : Vector{Int}[]
  out = Vector{Int}[]
  for first in 1:(n - j + 1), rest in compositions(n - first, j - 1)
    push!(out, [first; rest])
  end
  return out
end

@testset "peeling recursion equals composition sums" begin
  H = drive()
  N = 4
  vv = floquet_expansion(H, VanVleck(), N)

  function refH(n, j)
    acc = zero(H)
    for ks in compositions(n, j)
      term = H
      for k in Iterators.reverse(ks)
        term = comm(vv.K[k], term)
      end
      acc = acc + term
    end
    return (im^j * (1 // factorial(j))) * acc
  end

  function refKdot(n, j)
    acc = zero(H)
    for ks in compositions(n + 1, j + 1)
      term = vv.Kdot[ks[1]]
      for k in Iterators.reverse(ks[2:end])
        term = comm(vv.K[k], term)
      end
      acc = acc + term
    end
    return (im^j * (1 // factorial(j + 1))) * acc
  end

  for n in 0:(N - 1), j in 0:n
    node = FloquetExpansions.weightH(j) * vv.dressedH[FloquetExpansions.triindex(n, j)]
    @test agrees(node, refH(n, j))
    if j >= 1
      nodeK =
        FloquetExpansions.weightKdot(j) * vv.dressedKdot[FloquetExpansions.triindex(n, j)]
      @test agrees(nodeK, refKdot(n, j))
    end
  end
end

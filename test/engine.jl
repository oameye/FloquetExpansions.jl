using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra
using Symbolics: Symbolics

include(joinpath(@__DIR__, "helpers", "shared.jl"))

h = FockSpace(:cavity)
a = Destroy(h, :a)
@variables w::Real g::Real ξ::Real t::Real

function drive()
  return PeriodicGenerator(
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
function agrees(X::PeriodicGenerator, Y::PeriodicGenerator)
  return all(
    maxcoeff(SQA.simplify(X[l] - Y[l]), SUBS) < 1.0e-12 for l in union(keys(X), keys(Y))
  )
end

@testset "the spec's closed forms, orders 0 to 2" begin
  H = drive()
  Ms = [-2, -1, 1, 2]
  vv = floquet_expansion(H, VanVleck(), 3)

  @test vanishes(effective_component(vv, 0) - H[0])

  K1 = micromotion(vv, 1)
  @test sort!(collect(keys(K1))) == Ms
  @test all(vanishes(K1[m] - (im / (m * w)) * H[m]) for m in Ms)

  @test vanishes(
    effective_component(vv, 1) -
    sum((-1 // 2) * (1 / (m * w)) * comm(H[m], H[-m]) for m in Ms),
  )

  K2 = micromotion(vv, 2)
  function oracleK2(m)
    acc = comm(H[m], H[0]) * (1 / (m^2 * w^2))
    for mp in Ms
      mp == m && continue
      acc = acc + comm(H[mp], H[m - mp]) * (1 / (2 * m * mp * w^2))
    end
    return -im * acc
  end
  @test iszero(time_average(K2))
  @test all(vanishes(K2[m] - oracleK2(m)) for m in -4:4 if m != 0)
  @test_throws ArgumentError micromotion(vv, 0)
  @test_throws ArgumentError micromotion(vv, 3)

  oracleHeff2 = zero(SQA.QAdd)
  for m in Ms
    oracleHeff2 = oracleHeff2 + comm(comm(H[-m], H[0]), H[m]) * (1 / (2 * m^2 * w^2))
    for mp in Ms
      mp == m && continue
      oracleHeff2 =
        oracleHeff2 + comm(comm(H[-m], H[m - mp]), H[mp]) * (1 / (3 * m * mp * w^2))
    end
  end
  @test vanishes(effective_component(vv, 2) - oracleHeff2)
end

@testset "truncation follows the spec, X^[N] = sum_{k<N}" begin
  H = drive()

  rwa = floquet_expansion(H, VanVleck(), 1)
  @test vanishes(effective_generator(rwa) - H[0])
  @test iszero(micromotion(rwa))

  for N in 2:4
    lo = floquet_expansion(H, VanVleck(), N - 1)
    hi = floquet_expansion(H, VanVleck(), N)
    @test all(
      vanishes(effective_component(hi, n) - effective_component(lo, n)) for n in 0:(N - 2)
    )
    @test agrees(micromotion(hi) - micromotion(lo), micromotion(hi, N - 1))
  end
end

@testset "the kick round-trips through the collector" begin
  vv = floquet_expansion(drive(), VanVleck(), 3)
  K, back = micromotion(vv), harmonics(micromotion(vv)(t), w, t)
  @test sort!(collect(keys(back))) == sort!(collect(keys(K)))
  @test all(vanishes(back[l] - K[l]) for l in keys(K))
end

@testset "a non-Hermitian Hamiltonian drive is refused at ingest" begin
  bad = PeriodicGenerator(Dict(1 => a, -1 => a), w)
  @test_throws ArgumentError floquet_expansion(bad, VanVleck(), 2)
  @test_throws ArgumentError floquet_expansion(drive(), VanVleck(), 0)
end

@testset "the parsing entry point agrees with the harmonic one" begin
  H = 1 * (a' * a) + g * cos(w * t) * (a + a')
  viaparse = floquet_expansion(H, w, t, VanVleck(), 3)
  viaharm = floquet_expansion(harmonics(H, w, t), VanVleck(), 3)
  @test vanishes(effective_generator(viaparse) - effective_generator(viaharm))
end

function compositions(n::Int, j::Int)
  j == 0 && return n == 0 ? [Int[]] : Vector{Int}[]
  out = Vector{Int}[]
  for first in 1:(n - j + 1), rest in compositions(n - first, j - 1)
    push!(out, [first; rest])
  end
  return out
end

function reference_expansion(H::PeriodicGenerator, wd::Symbolics.Num, N::Int)
  kicks = typeof(H)[]
  kick_derivatives = typeof(H)[]
  effective = typeof(time_average(H))[]

  for n in 0:(N - 1)
    resolvent = zero(H)

    for j in 0:n
      contribution = zero(H)
      for ks in compositions(n, j)
        term = H
        for k in Iterators.reverse(ks)
          term = comm(kicks[k], term)
        end
        contribution = contribution + term
      end
      resolvent = resolvent + (im^j * (1 // factorial(j))) * contribution
    end

    for j in 1:n
      contribution = zero(H)
      for ks in compositions(n + 1, j + 1)
        term = kick_derivatives[ks[1]]
        for k in Iterators.reverse(ks[2:end])
          term = comm(kicks[k], term)
        end
        contribution = contribution + term
      end
      resolvent = resolvent - (im^j * (1 // factorial(j + 1))) * contribution
    end

    resolvent = SQA.simplify(resolvent)
    push!(effective, time_average(resolvent))

    if n < N - 1
      average = time_average(resolvent)
      oscillatory = resolvent - PeriodicGenerator(Dict(0 => average), wd)
      next_kick = SQA.simplify(antiderivative(oscillatory, VanVleck()))
      push!(kicks, next_kick)
      push!(kick_derivatives, derivative(next_kick))
    end
  end

  return effective, kicks
end

@testset "expansion matches composition oracle" begin
  H = drive()
  N = 4
  vv = floquet_expansion(H, VanVleck(), N)
  expected_effective, expected_kicks = reference_expansion(H, w, N)

  for n in 0:(N - 1)
    expected = w^(-n) * expected_effective[n + 1]
    @test maxcoeff(effective_component(vv, n) - expected, SUBS) < 1.0e-12
  end

  for n in 1:(N - 1)
    @test agrees(micromotion(vv, n), w^(-n) * expected_kicks[n])
  end
end

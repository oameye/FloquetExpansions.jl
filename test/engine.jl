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

@testset "the spec's closed forms, orders 0 to 2" begin
  H = drive()
  Ms = [-2, -1, 1, 2]
  vv = floquet_expansion(H, VanVleck(), 3)

  @test vanishes(effective_hamiltonian(vv, 0) - H[0])

  K1 = kick_operator(vv, 1)
  @test sort!(collect(keys(K1))) == Ms
  @test all(vanishes(K1[m] - (im / (m * w)) * H[m]) for m in Ms)

  @test vanishes(
    effective_hamiltonian(vv, 1) -
    sum((-1 // 2) * (1 / (m * w)) * comm(H[m], H[-m]) for m in Ms),
  )

  K2 = kick_operator(vv, 2)
  function oracleK2(m)
    acc = comm(H[m], H[0]) * (1 / (m^2 * w^2))
    for mp in Ms
      mp == m && continue
      acc = acc + comm(H[mp], H[m - mp]) * (1 / (2 * m * mp * w^2))
    end
    return -im * acc
  end
  @test !haskey(K2.components, 0)
  @test all(vanishes(K2[m] - oracleK2(m)) for m in -4:4 if m != 0)
  @test_throws ArgumentError kick_operator(vv, 0)
  @test_throws ArgumentError kick_operator(vv, 3)

  oracleHeff2 = zero(SQA.QAdd)
  for m in Ms
    oracleHeff2 = oracleHeff2 + comm(comm(H[-m], H[0]), H[m]) * (1 / (2 * m^2 * w^2))
    for mp in Ms
      mp == m && continue
      oracleHeff2 =
        oracleHeff2 + comm(comm(H[-m], H[m - mp]), H[mp]) * (1 / (3 * m * mp * w^2))
    end
  end
  @test vanishes(effective_hamiltonian(vv, 2) - oracleHeff2)
end

@testset "truncation follows the spec, X^[N] = sum_{k<N}" begin
  H = drive()

  rwa = floquet_expansion(H, VanVleck(), 1)
  @test vanishes(effective_hamiltonian(rwa) - H[0])
  @test iszero(kick_operator(rwa))

  for N in 2:4
    lo = floquet_expansion(H, VanVleck(), N - 1)
    hi = floquet_expansion(H, VanVleck(), N)
    @test all(
      vanishes(effective_hamiltonian(hi, n) - effective_hamiltonian(lo, n)) for
      n in 0:(N - 2)
    )
    @test agrees(kick_operator(hi) - kick_operator(lo), kick_operator(hi, N - 1))
  end
end

@testset "the kick round-trips through the collector" begin
  vv = floquet_expansion(drive(), VanVleck(), 3)
  K, back = kick_operator(vv), harmonics(kick_operator(vv, t), w, t)
  @test sort!(collect(keys(back))) == sort!(collect(keys(K)))
  @test all(vanishes(back[l] - K[l]) for l in keys(K))
end

@testset "a non-Hermitian drive is refused at ingest" begin
  bad = PeriodicOperator(Dict(1 => a, -1 => a), w)
  @test_throws ArgumentError floquet_expansion(bad, VanVleck(), 2)
  @test_throws ArgumentError floquet_expansion(drive(), VanVleck(), 0)
end

@testset "the parsing entry point agrees with the harmonic one" begin
  H = 1 * (a' * a) + g * cos(w * t) * (a + a')
  viaparse = floquet_expansion(H, w, t, VanVleck(), 3)
  viaharm = floquet_expansion(harmonics(H, w, t), VanVleck(), 3)
  @test vanishes(effective_hamiltonian(viaparse) - effective_hamiltonian(viaharm))
end

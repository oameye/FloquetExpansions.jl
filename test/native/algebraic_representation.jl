using Test
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions
const Q = Complex{Rational{BigInt}}
const R = Rational{BigInt}

@testset "cutoff-free operator algebra obeys the canonical relations" begin
  algebra = FE.OperatorAlgebra{Q}([
    FE.boson_site(R), FE.phase_site(R), FE.level_site(R[1, 2])
  ])
  a = FE.algebra_operator(algebra, [([0, 1, 0, 0, 0, 0], 1)])
  x = FE.algebra_operator(algebra, [([0, 0, 1, 0, 0, 0], 1)])
  p = FE.algebra_operator(algebra, [([0, 0, 0, 1, 0, 0], 1)])
  σ = FE.algebra_operator(algebra, [([0, 0, 0, 0, 1, 2], 1)])
  identity = FE.algebra_identity(algebra)
  @test a * adjoint(a) - adjoint(a) * a == identity
  @test x * p - p * x == im * identity
  @test adjoint(x) == x
  @test adjoint(p) == p
  # With Hilbert metric (1, 2) the adjoint of |1)(2| is (g_1 / g_2) |2)(1|.
  @test adjoint(σ) == FE.algebra_operator(algebra, [([0, 0, 0, 0, 2, 1], 1 // 2)])
  @test σ * adjoint(σ) + adjoint(σ) * σ == (1 // 2) * identity
  @test a * x == x * a
end

function cutoff_free_kerr(N, algorithm)
  algebra = FE.OperatorAlgebra{Q}([FE.boson_site(R)])
  rep = FE.AlgebraicLiouvilleRepresentation(algebra, 4, 1)
  a = FE.algebra_operator(algebra, [([0, 1], 1)])
  n = adjoint(a) * a
  x = a + adjoint(a)
  c = FE.frame_coordinates(rep, a)
  Z = zeros(Q, length(rep.frame), length(rep.frame))
  F = 2 // 5 + im // 5
  L = Dict(
    0 => FE.native_gksl(rep, (1 // 2) * n + (3 // 10) * n * n, (4 // 5) * c * c'),
    1 => FE.native_gksl(rep, F * x, Z),
    -1 => FE.native_gksl(rep, conj(F) * x, Z),
  )
  leading = reshape(c, :, 1)
  inverse = FE.AlgebraicChargeInverse(rep, L[0], leading, 6)
  return FE.native_recurrence(algorithm, rep, inverse, L, leading, Q[4 // 5], N, 1e-8), rep
end

@testset "cutoff-free driven Kerr: BF equals HD with no truncation births through order 2" begin
  bf, rep = cutoff_free_kerr(2, FE.BlochFeshbach())
  hd, _ = cutoff_free_kerr(2, FE.HoriDeprit())
  @test bf.E == hd.E
  @test bf.H == hd.H
  @test length(bf.channels) == 1
  @test all(iszero, bf.S)
  for order in 0:2
    @test FE.native_kossakowski(rep, bf.E[order + 1]) ==
      FE.gram_coefficient(bf.channels, order, length(rep.frame))
  end
end

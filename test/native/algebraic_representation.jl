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
  @test adjoint(σ) == FE.algebra_operator(algebra, [([0, 0, 0, 0, 2, 1], 1 // 2)])
  @test σ * adjoint(σ) + adjoint(σ) * σ == (1 // 2) * identity
  @test a * x == x * a
end

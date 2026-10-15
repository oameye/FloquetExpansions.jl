using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

fock = FockSpace(:gksl_floquet)
a = Destroy(fock, :a)
@variables Δ::Real ω::Real t::Real

@testset "Floquet GKSL coordinate accessors" begin
  H = Δ * a' * a
  L = SQA.SuperOperator(hamiltonian_action(H) + dissipator(a + im * a^2))
  frame = JumpBasis(a, a^2)
  expansion = floquet_expansion(L, ω, t, VanVleck(), 1)

  d = @inferred kossakowski(expansion, frame)
  d0 = @inferred kossakowski_component(expansion, frame, 0)
  H_eff = @inferred hamiltonian(expansion)
  H0 = @inferred hamiltonian_component(expansion, 0)

  @test d == d0
  @test d[1, 2] == -im
  @test iszero(SQA.simplify(H_eff - H))
  @test iszero(SQA.simplify(H0 - H))
end

@testset "nonzero first-order Kossakowski component" begin
  coherent_harmonic = SQA.SuperOperator(hamiltonian_action(a'^2 * a^2))
  dissipative_harmonic = SQA.SuperOperator(dissipator(a))
  quadrature = im * (1 // 2) * dissipative_harmonic
  generator = PeriodicGenerator(
    Dict(1 => coherent_harmonic - quadrature, -1 => coherent_harmonic + quadrature), ω
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)
  frame = JumpBasis(a, a' * a^2)

  d1 = @inferred kossakowski_component(expansion, frame, 1)

  @test all(
    iszero(SQA.simplify(d1[i, j] - conj(d1[j, i]))) for i in axes(d1, 1), j in axes(d1, 2)
  )
  @test iszero(SQA.simplify(d1[1, 1]))
  @test iszero(SQA.simplify(d1[2, 2]))
  @test iszero(SQA.simplify(d1[1, 2] + 2im / ω))
  @test iszero(SQA.simplify(d1[2, 1] - 2im / ω))
end

@testset "Hamiltonian expansions retain their coherent accessors" begin
  H = Δ * a' * a
  expansion = floquet_expansion(harmonics(H, ω, t), VanVleck(), 1)

  @test @inferred(hamiltonian(expansion)) == effective_generator(expansion)
  @test @inferred(hamiltonian_component(expansion, 0)) == effective_component(expansion, 0)
end

using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

fock = FockSpace(:gksl_coordinates)
a = Destroy(fock, :a)
@variables Δ::Real ω::Real t::Real

@testset "DissipativeFrame is ordered and independent modulo identity" begin
  frame = @inferred DissipativeFrame(a, a^2)

  @test frame == DissipativeFrame((a, a^2))
  @test frame != DissipativeFrame(a^2, a)
  @test_throws ArgumentError DissipativeFrame(a, 2a)
  @test_throws ArgumentError DissipativeFrame(one(a), a)
end

@testset "exact complex Kossakowski extraction" begin
  H = Δ * a' * a
  frame = DissipativeFrame(a, a^2)
  L = liouvillian(H; channels=(collapse(a + im * a^2),))

  d = @inferred kossakowski(L, frame)
  extracted_H = @inferred hamiltonian(L, frame)

  @test size(d) == (2, 2)
  @test d[1, 1] == 1
  @test d[1, 2] == -im
  @test d[2, 1] == im
  @test d[2, 2] == 1
  @test all(d[i, j] == conj(d[j, i]) for i in axes(d, 1), j in axes(d, 2))
  @test iszero(SQA.simplify(extracted_H - H))
  @test iszero(SQA.simplify(hamiltonian(L) - H))
end

@testset "nonorthogonal rational dissipative frame" begin
  f1 = a + a^2
  f2 = a - a^2
  frame = DissipativeFrame(f1, f2)
  H = (3 // 5) * a' * a
  L = liouvillian(H; channels=(collapse(2 * f1 + 3 * f2),))

  d = kossakowski(L, frame)

  @test d[1, 1] == 4
  @test d[1, 2] == 6
  @test d[2, 1] == 6
  @test d[2, 2] == 9
  @test hamiltonian_action(hamiltonian(L, frame)) == hamiltonian_action(H)
end

@testset "Hamiltonian is recovered modulo the identity" begin
  frame = DissipativeFrame(a)
  H = a' * a + 7 * one(a)
  L = liouvillian(H; channels=(collapse(a),))

  extracted_H = hamiltonian(L, frame)

  @test iszero(SQA.simplify(extracted_H - a' * a))
  @test kossakowski(L, frame)[1, 1] == 1
end

@testset "NLevelSpace completeness is expanded before coordinates" begin
  finite = NLevelSpace(:gksl_finite, 2)
  sigma11 = Transition(finite, :sigma, 1, 1)
  frame = DissipativeFrame(sigma11)
  L = liouvillian(zero(SQA.QAdd); channels=(collapse(sigma11),))

  d = @inferred kossakowski(L, frame)

  @test d[1, 1] == 1
  @test iszero(hamiltonian(L, frame))
end

@testset "bosonic monomial frame needs no Hilbert cutoff" begin
  frame = DissipativeFrame(a, a^2, a' * a^2)
  channel = a + 2a^2 - im * a' * a^2
  L = liouvillian(zero(SQA.QAdd); channels=(collapse(channel),))

  d = kossakowski(L, frame)

  @test size(d) == (3, 3)
  @test d[1, 1] == 1
  @test d[1, 2] == 2
  @test d[1, 3] == im
  @test d[2, 1] == 2
  @test d[2, 2] == 4
  @test d[2, 3] == 2im
  @test d[3, 1] == -im
  @test d[3, 2] == -2im
  @test d[3, 3] == 1
end

@testset "incomplete frame is rejected" begin
  L = liouvillian(zero(SQA.QAdd); channels=(collapse(a + a^2),))
  frame = DissipativeFrame(a)

  @test_throws ArgumentError kossakowski(L, frame)
  @test_throws ArgumentError hamiltonian(L, frame)
end

@testset "Floquet GKSL coordinate accessors" begin
  H = Δ * a' * a
  L = liouvillian(H; channels=(collapse(a + im * a^2),))
  frame = DissipativeFrame(a, a^2)
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

@testset "Hamiltonian expansions retain their coherent accessors" begin
  H = Δ * a' * a
  expansion = floquet_expansion(harmonics(H, ω, t), VanVleck(), 1)

  @test @inferred(hamiltonian(expansion)) == effective_generator(expansion)
  @test @inferred(hamiltonian_component(expansion, 0)) == effective_generator(expansion, 0)
end

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
  shifted = DissipativeFrame(a + 3 * one(a), a^2)

  @test frame == DissipativeFrame((a, a^2))
  @test frame == DissipativeFrame([a, a^2])
  @test frame == shifted
  @test isequal(frame, shifted)
  @test hash(frame) == hash(shifted)
  @test frame != DissipativeFrame(a^2, a)
  @test_throws ArgumentError DissipativeFrame()
  @test_throws ArgumentError DissipativeFrame(a, 2a)
  @test_throws ArgumentError DissipativeFrame(one(a), a)
  @test_throws ArgumentError DissipativeFrame(Any[a, 1])
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

@testset "nonorthogonal frames transform Kossakowski coordinates by congruence" begin
  f1 = a + a^2
  f2 = a - a^2
  rotated_frame = DissipativeFrame(f1, f2)
  native_frame = DissipativeFrame(a, a^2)
  H = (3 // 5) * a' * a
  L = liouvillian(H; channels=(collapse(2 * f1 + 3 * f2),))

  rotated = kossakowski(L, rotated_frame)
  native = kossakowski(L, native_frame)
  transform = [1 1; 1 -1]

  @test rotated[1, 1] == 4
  @test rotated[1, 2] == 6
  @test rotated[2, 1] == 6
  @test rotated[2, 2] == 9
  @test native == transform * rotated * transform'
  @test hamiltonian_action(hamiltonian(L, rotated_frame)) == hamiltonian_action(H)
end

@testset "Hamiltonian is recovered modulo the identity" begin
  frame = DissipativeFrame(a)
  H = a' * a + 7 * one(a)
  L = liouvillian(H; channels=(collapse(a),))

  extracted_H = hamiltonian(L, frame)

  @test iszero(SQA.simplify(extracted_H - a' * a))
  @test kossakowski(L, frame)[1, 1] == 1
end

@testset "pure Hamiltonian Liouvillian has a zero dissipative sector" begin
  H = Δ * a' * a + 5 * one(a)
  L = liouvillian(H)
  frame = DissipativeFrame(a, a^2)

  extracted_H = @inferred hamiltonian(L)
  d = @inferred kossakowski(L, frame)

  @test iszero(SQA.simplify(extracted_H - Δ * a' * a))
  @test all(iszero, d)
end

@testset "inactive frame directions remain exactly dark" begin
  frame = DissipativeFrame(a, a^2)
  L = liouvillian(0 * a; channels=(collapse(a),))

  d = kossakowski(L, frame)

  @test d[1, 1] == 1
  @test iszero(d[1, 2])
  @test iszero(d[2, 1])
  @test iszero(d[2, 2])
end

@testset "adjoint thermal channels retain independent rates" begin
  frame = DissipativeFrame(a, a')
  L = liouvillian(0 * a; channels=(collapse(a), jump(a', 2)))

  d = kossakowski(L, frame)

  @test d[1, 1] == 1
  @test iszero(d[1, 2])
  @test iszero(d[2, 1])
  @test d[2, 2] == 2
  @test iszero(hamiltonian(L, frame))
end

@testset "identity shifts move only into the Hamiltonian gauge" begin
  shifted_jump = a + 3 * one(a)
  frame = DissipativeFrame(shifted_jump)
  L = liouvillian(0 * a; channels=(collapse(shifted_jump),))
  expected_H = (3 // 2) * im * (a - a')

  d = kossakowski(L, frame)
  extracted_H = hamiltonian(L, frame)

  @test d[1, 1] == 1
  @test hamiltonian_action(extracted_H) == hamiltonian_action(expected_H)
end

@testset "NLevelSpace completeness is expanded before coordinates" begin
  finite = NLevelSpace(:gksl_finite, 2)
  sigma11 = Transition(finite, :sigma, 1, 1)
  frame = DissipativeFrame(sigma11)
  L = liouvillian(0 * sigma11; channels=(collapse(sigma11),))

  d = @inferred kossakowski(L, frame)

  @test d[1, 1] == 1
  @test iszero(hamiltonian(L, frame))
end

@testset "Pauli algebra closes inside exact GKSL coordinates" begin
  pauli = PauliSpace(:gksl_pauli)
  sigma_x = Pauli(pauli, :sigma, 1)
  sigma_y = Pauli(pauli, :sigma, 2)
  sigma_z = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(sigma_x, sigma_y)
  H = Δ * sigma_z
  L = liouvillian(H; channels=(collapse(sigma_x + im * sigma_y),))

  d = kossakowski(L, frame)

  @test d[1, 1] == 1
  @test d[1, 2] == -im
  @test d[2, 1] == im
  @test d[2, 2] == 1
  @test hamiltonian_action(hamiltonian(L, frame)) == hamiltonian_action(H)
end

@testset "canonical phase-space algebra preserves dissipative coordinates" begin
  phase = PhaseSpace(:gksl_phase)
  x = Position(phase, :x)
  p = Momentum(phase, :p)
  frame = DissipativeFrame(x, p)
  L = liouvillian(0 * x; channels=(collapse(x + im * p),))

  d = kossakowski(L, frame)

  @test d[1, 1] == 1
  @test d[1, 2] == -im
  @test d[2, 1] == im
  @test d[2, 2] == 1
  @test iszero(hamiltonian(L, frame))
end

@testset "mixed product-space directions remain independent" begin
  mixed = FockSpace(:gksl_mixed_boson) ⊗ PauliSpace(:gksl_mixed_spin)
  b = Destroy(mixed, :b)
  sigma_x = Pauli(mixed, :sigma, 1)
  frame = DissipativeFrame(b, sigma_x)
  L = liouvillian(0 * b; channels=(collapse(b + 2 * sigma_x),))

  d = kossakowski(L, frame)

  @test d[1, 1] == 1
  @test d[1, 2] == 2
  @test d[2, 1] == 2
  @test d[2, 2] == 4
  @test iszero(hamiltonian(L, frame))
end

@testset "bosonic monomial frame needs no Hilbert cutoff" begin
  frame = DissipativeFrame(a, a^2, a' * a^2)
  channel = a + 2a^2 - im * a' * a^2
  L = liouvillian(0 * a; channels=(collapse(channel),))

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
  L = liouvillian(0 * a; channels=(collapse(a + a^2),))
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

@testset "nonzero first-order Kossakowski component" begin
  coherent_harmonic = hamiltonian_action(a'^2 * a^2)
  dissipative_harmonic = dissipator(a)
  generator = PeriodicGenerator(Dict(-1 => coherent_harmonic, 1 => dissipative_harmonic), ω)
  expansion = floquet_expansion(generator, VanVleck(), 2)
  frame = DissipativeFrame(a, a' * a^2)

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
  @test @inferred(hamiltonian_component(expansion, 0)) == effective_generator(expansion, 0)
end

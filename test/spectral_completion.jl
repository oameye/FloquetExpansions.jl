using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

function spectral_matrix_equal(left, right)
  size(left) == size(right) || return false
  return all(iszero(SQA.simplify(left[index] - right[index])) for index in eachindex(left))
end

fock = FockSpace(:spectral_completion)
a = Destroy(fock, :a)
@variables ω::Real t::Real γ::Real Ω::Real

@testset "spectral completion preserves a diagonal rate branch" begin
  frame = DissipativeFrame(a)
  expansion = floquet_expansion(0 * a, ω, t, VanVleck(), 1; channels=(jump(a, γ),))
  completion = @inferred positive_completion(expansion, Spectral(), frame)
  spectral = factorization(completion)

  @test spectral isa SpectralFactorization
  @test spectral.onsets == [0]
  @test spectral.puiseux == [false]
  @test dissipative_frame(completion) == frame
  @test spectral_matrix_equal(kossakowski(completion), kossakowski(expansion, frame))
  @test effective_component(completion, 0) == effective_component(expansion, 0)
  @test micromotion(completion) == micromotion(expansion)
  @test liouvillian(hamiltonian(completion); channels=channels(completion)) ==
    effective_generator(completion)
end

@testset "adapted-frame driven-qubit spectral completion" begin
  pauli = PauliSpace(:spectral_pauli)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  expansion = floquet_expansion(
    Ω * cos(ω * t) * σx, ω, t, VanVleck(), 3; channels=(collapse(σz),)
  )
  completion = @inferred positive_completion(expansion, Spectral(), frame)
  spectral = factorization(completion)

  @test spectral isa SpectralFactorization
  @test length(spectral.rates) == 2
  @test length(spectral.vectors) == 2
  @test all(length(vector) == 2 for vector in spectral.vectors)
  for n in 0:2
    @test effective_component(completion, n) == effective_component(expansion, n)
    @test spectral_matrix_equal(
      kossakowski_component(completion, n), kossakowski_component(expansion, frame, n)
    )
  end
  @test micromotion(completion) == micromotion(expansion)
  @test liouvillian(hamiltonian(completion); channels=channels(completion)) ==
    effective_generator(completion)
end

@testset "distinct leading rates generate perturbative branch mixing" begin
  pauli = PauliSpace(:spectral_mixing)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  coherent_rotation = hamiltonian_action(σx)
  diagonal_difference = dissipator(σy) - dissipator(σz)
  leading = dissipator(σz) + 4 * dissipator(σy)
  generator = PeriodicGenerator(
    Dict(0 => leading, 1 => coherent_rotation, -1 => diagonal_difference), ω
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)
  mixing = kossakowski_component(expansion, frame, 1)

  @test !iszero(SQA.simplify(mixing[1, 2]))
  @test !iszero(SQA.simplify(mixing[2, 1]))

  completion = positive_completion(expansion, Spectral(), frame)
  spectral = factorization(completion)
  @test spectral.onsets == [0, 0]
  @test any(
    !iszero(SQA.simplify(spectral.vectors[column][row]))
    for (column, row) in ((1, 2), (2, 1))
  )
  for n in 0:1
    @test spectral_matrix_equal(
      kossakowski_component(completion, n), kossakowski_component(expansion, frame, n)
    )
  end
end

@testset "degenerate leading sector with retained mixing is rejected" begin
  pauli = PauliSpace(:spectral_degenerate)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  coherent_rotation = hamiltonian_action(σx)
  diagonal_difference = dissipator(σy) - dissipator(σz)
  leading = dissipator(σz) + dissipator(σy)
  generator = PeriodicGenerator(
    Dict(0 => leading, 1 => coherent_rotation, -1 => diagonal_difference), ω
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)
  mixing = kossakowski_component(expansion, frame, 1)

  @test !iszero(SQA.simplify(mixing[1, 2]))
  @test !iszero(SQA.simplify(mixing[2, 1]))
  @test_throws ArgumentError positive_completion(expansion, Spectral(), frame)
end

@testset "odd positive rate onset remains a valid spectral channel" begin
  pauli = PauliSpace(:spectral_odd_onset)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)

  coherent_rotation = hamiltonian_action(σx)
  cross_dissipator = dissipator(σy + σz) - dissipator(σy) - dissipator(σz)
  generator = PeriodicGenerator(
    Dict(0 => dissipator(σz), 1 => coherent_rotation, -1 => cross_dissipator), ω
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)
  completion = positive_completion(expansion, Spectral(), frame)
  spectral = factorization(completion)

  @test 1 in spectral.onsets
  odd_branch = findfirst(==(1), spectral.onsets)
  @test odd_branch !== nothing
  @test spectral.puiseux[odd_branch]
  @test liouvillian(hamiltonian(completion); channels=channels(completion)) ==
    effective_generator(completion)
  for n in 0:1
    @test spectral_matrix_equal(
      kossakowski_component(completion, n), kossakowski_component(expansion, frame, n)
    )
  end
end

@testset "Gram and spectral continuations agree through retained qubit order" begin
  pauli = PauliSpace(:spectral_gram_comparison)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  expansion = floquet_expansion(
    Ω * cos(ω * t) * σx, ω, t, VanVleck(), 3; channels=(collapse(σz),)
  )
  gram = positive_completion(expansion, Gram(), frame)
  spectral = positive_completion(expansion, Spectral(), frame)

  for n in 0:2
    retained = kossakowski_component(expansion, frame, n)
    @test spectral_matrix_equal(kossakowski_component(gram, n), retained)
    @test spectral_matrix_equal(kossakowski_component(spectral, n), retained)
  end
  @test effective_component(gram, 2) == effective_component(spectral, 2)
  @test micromotion(gram) == micromotion(spectral)
  @test liouvillian(hamiltonian(gram); channels=channels(gram)) == effective_generator(gram)
  @test liouvillian(hamiltonian(spectral); channels=channels(spectral)) ==
    effective_generator(spectral)
end

@testset "non-diagonal leading spectral frame is rejected" begin
  frame = DissipativeFrame(a, a^2)
  expansion = floquet_expansion(
    0 * a, ω, t, VanVleck(), 1; channels=(collapse(a + a^2), collapse(a + im * a^2))
  )
  @test_throws ArgumentError positive_completion(expansion, Spectral(), frame)
end

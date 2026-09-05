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

@testset "non-diagonal leading spectral frame is rejected" begin
  frame = DissipativeFrame(a, a^2)
  expansion = floquet_expansion(
    0 * a, ω, t, VanVleck(), 1; channels=(collapse(a + a^2), collapse(a + im * a^2))
  )
  @test_throws ArgumentError positive_completion(expansion, Spectral(), frame)
end

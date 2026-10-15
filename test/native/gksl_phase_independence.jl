using Test
using Symbolics: @variables
using FloquetExpansions
using FloquetExpansions: FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra

const FE = FloquetExpansions
const SQA = SecondQuantizedAlgebra

@variables ω::Real t::Real

function same(A::SQA.SuperOperator, B::SQA.SuperOperator)
  return iszero(SQA.simplify(SQA.expand_completeness(A - B)))
end
same(A, B) = iszero(SQA.simplify(FE.qadd(A) - FE.qadd(B)))

channel_map(entry) = SQA.SuperOperator(last(entry) * dissipator(SQA.jump(first(entry))))

phased_cosine() = ((3 + 4im) // 10) * expim(ω * t) + ((3 - 4im) // 10) * expim(-ω * t)
plain_cosine() = (1 // 2) * expim(ω * t) + (1 // 2) * expim(-ω * t)

function driven_kerr(drive)
  h = FockSpace(:gksl_phase)
  a = Destroy(h, :a)
  H = (1 // 2) * a' * a + (3 // 10) * a' * a * a' * a + (2 // 5) * (a + a') * drive
  return H, a
end

function dissipative_part(expansion)
  return sum(channel_map, channels(expansion); init=zero(SQA.SuperOperator))
end

function has_channel(expansion, expected::Lindbladian)
  return count(
    entry -> same(channel_map(entry), SQA.SuperOperator(expected)), channels(expansion)
  ) == 1
end

@testset "GKSLNormalForm is independent of the drive phase" begin
  H0, a = driven_kerr(plain_cosine())
  Hφ, _ = driven_kerr(phased_cosine())
  for algorithm in (HoriDeprit(), BlochFeshbach())
    gauge = GKSLNormalForm(; algorithm)
    loss = (4 // 5) * dissipator(a)
    plain = floquet_expansion(hamiltonian_action(H0) + loss, ω, t, gauge, 5)
    phased = floquet_expansion(hamiltonian_action(Hφ) + loss, ω, t, gauge, 5)
    @test length(channels(plain)) == length(channels(phased)) == 2
    @test same(hamiltonian(plain), hamiltonian(phased))
    @test same(effective_generator(plain), effective_generator(phased))
    @test same(dissipative_part(plain), dissipative_part(phased))
    for n in 0:4
      @test same(effective_component(plain, n), effective_component(phased, n))
    end
    induced = (72 // 3125) / ω^4 * dissipator(a * a - 2 * a' * a)
    @test has_channel(plain, induced)
    @test has_channel(phased, induced)
  end
end

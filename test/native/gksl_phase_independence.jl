using Test
using Symbolics: @variables
using FloquetExpansions
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions

@variables ω::Real t::Real

function same(A::Liouvillian, B::Liouvillian)
  return FE.liouvillian_iszero(FE.canonical_liouvillian(FE.SQA.simplify(A - B)))
end
same(A, B) = iszero(FE.SQA.simplify(FE.qadd(A) - FE.qadd(B)))

# cos(ωt + φ) with exp(iφ) = (3 + 4i) / 5, written with exact Gaussian-rational harmonics
phased_cosine() = ((3 + 4im) // 10) * expim(ω * t) + ((3 - 4im) // 10) * expim(-ω * t)
plain_cosine() = (1 // 2) * expim(ω * t) + (1 // 2) * expim(-ω * t)

function driven_kerr(drive)
  h = FockSpace(:gksl_phase)
  a = Destroy(h, :a)
  H = (1 // 2) * a' * a + (3 // 10) * a' * a * a' * a + (2 // 5) * (a + a') * drive
  return H, a
end

function dissipative_part(expansion)
  return liouvillian(0 * hamiltonian(expansion); channels=Tuple(channels(expansion)))
end

@testset "GKSLNormalForm is independent of the drive phase" begin
  H0, a = driven_kerr(plain_cosine())
  Hφ, _ = driven_kerr(phased_cosine())
  for algorithm in (HoriDeprit(), BlochFeshbach())
    gauge = GKSLNormalForm(; algorithm)
    plain = floquet_expansion(H0, ω, t, gauge, 5; channels=(jump(a, 4 // 5),))
    phased = floquet_expansion(Hφ, ω, t, gauge, 5; channels=(jump(a, 4 // 5),))
    @test length(channels(plain)) == length(channels(phased)) == 2
    @test same(hamiltonian(plain), hamiltonian(phased))
    @test same(effective_generator(plain), effective_generator(phased))
    @test same(dissipative_part(plain), dissipative_part(phased))
    for n in 0:4
      @test same(effective_component(plain, n), effective_component(phased, n))
    end
    plain_rates = [FE.SQA.simplify(c.rate * ω^4) for c in channels(plain)][2]
    phased_rates = [FE.SQA.simplify(c.rate * ω^4) for c in channels(phased)][2]
    @test isequal(plain_rates, 72 // 3125)
    @test isequal(phased_rates, 72 // 3125)
    @test isequal(channels(plain)[1].rate, channels(phased)[1].rate)
  end
end

using Test
using Symbolics: @variables
using FloquetExpansions
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions

@variables ω::Real t::Real Δ::Real χ::Real κ::Real F::Real

const VALUES = Dict(Δ => 1 // 2, χ => 3 // 10, κ => 4 // 5, F => 2 // 5)

numeric(c) = ComplexF64(FE.SQA.to_complex(FE.SQA.substitute(c, Dict(VALUES..., ω => 10))))

function evaluated(L::Liouvillian)
  result = zero(Liouvillian)
  for (left, right, c) in FE.terms(L)
    result = result + numeric(c) * FE.action(left, right)
  end
  return result
end

function nearly_same(A::Liouvillian, B::Liouvillian)
  difference = FE.canonical_liouvillian(FE.SQA.simplify(evaluated(A) - evaluated(B)))
  return all(abs(numeric(c)) <= 1e-12 for (_, _, c) in FE.terms(difference))
end

function kerr(drive)
  h = FockSpace(:c)
  a = Destroy(h, :a)
  return Δ * a' * a + χ * a' * a * a' * a + drive * (a + a') * cos(ω * t), a
end

@testset "symbolic GKSLNormalForm without a static gauge is the Gram completion" begin
  H, a = kerr(F)
  symbolic = floquet_expansion(H, ω, t, GKSLNormalForm(), 4; channels=(jump(a, κ),))
  @test getfield(symbolic, :completion) isa FE.PositiveCompletion
  @test any(c -> isequal(FE.SQA.simplify(c - κ), 0), positivity_conditions(symbolic))
  @test_throws ArgumentError positive_completion(symbolic, Gram())

  Hn = 1 // 2 * a' * a + 3 // 10 * a' * a * a' * a + 2 // 5 * (a + a') * cos(ω * t)
  exact = floquet_expansion(Hn, ω, t, GKSLNormalForm(), 4; channels=(jump(a, 4 // 5),))
  @test nearly_same(effective_generator(symbolic), effective_generator(exact))
  for n in 0:3
    @test nearly_same(effective_component(symbolic, n), effective_component(exact, n))
  end
end

@testset "symbolic GKSLNormalForm reports an order that needs a static gauge" begin
  H, a = kerr(F)
  error = try
    floquet_expansion(H, ω, t, GKSLNormalForm(), 5; channels=(jump(a, κ),))
  catch caught
    caught
  end
  @test error isa ArgumentError
  @test occursin("static gauge at order 4", sprint(showerror, error))
end

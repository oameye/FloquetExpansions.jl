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

function nearly_same(A::FE.SQA.QAdd, B::FE.SQA.QAdd)
  return nearly_same(FE.action(A, one(A)), FE.action(B, one(B)))
end

function nearly_same_micromotion(A::FloquetExpansion, B::FloquetExpansion, n::Int)
  m, r = micromotion(A, n), micromotion(B, n)
  return all(
    nearly_same(
      get(m.components, k, zero(m.zero_component)),
      get(r.components, k, zero(r.zero_component)),
    ) for k in union(keys(m.components), keys(r.components))
  )
end

symbolically_same(x, y) = isequal(FE.SQA.simplify(x - y), 0)

function same_liouvillian(A::Liouvillian, B::Liouvillian)
  return FE.liouvillian_iszero(FE.canonical_liouvillian(FE.SQA.simplify(A - B)))
end

function newborn_channel(expansion::FloquetExpansion)
  return only(c for c in channels(expansion) if !symbolically_same(c.rate, κ))
end

@testset "symbolic GKSLNormalForm reconstructs the order-four static gauge" begin
  H, a = kerr(F)
  symbolic = floquet_expansion(H, ω, t, GKSLNormalForm(), 5; channels=(jump(a, κ),))
  @test getfield(symbolic, :completion) isa FE.NativeRealization
  @test length(channels(symbolic)) == 2
  newborn = newborn_channel(symbolic)
  @test symbolically_same(newborn.rate, convert(FE.SQA.CNum, 2 * χ^2 * κ * F^2 / ω^4))
  @test iszero(FE.SQA.simplify(newborn.operator - (a * a - 2 * a' * a)))
  conditions = positivity_conditions(symbolic)
  for c in (Δ, χ, κ, F, 2 * χ^2 * κ * F^2)
    @test any(d -> symbolically_same(d, convert(FE.SQA.CNum, c)), conditions)
  end

  Hn = 1 // 2 * a' * a + 3 // 10 * a' * a * a' * a + 2 // 5 * (a + a') * cos(ω * t)
  exact = floquet_expansion(Hn, ω, t, GKSLNormalForm(), 5; channels=(jump(a, 4 // 5),))
  @test nearly_same(effective_generator(symbolic), effective_generator(exact))
  @test nearly_same(hamiltonian(symbolic), hamiltonian(exact))
  for n in 0:4
    @test nearly_same(effective_component(symbolic, n), effective_component(exact, n))
  end
  for n in 1:4
    @test nearly_same_micromotion(symbolic, exact, n)
  end
  @test sort(numeric.(c.rate for c in channels(symbolic)); by=real) ≈
    sort(numeric.(c.rate for c in channels(exact)); by=real)

  bloch = floquet_expansion(
    H, ω, t, GKSLNormalForm(; algorithm=BlochFeshbach()), 5; channels=(jump(a, κ),)
  )
  for n in 0:4
    @test same_liouvillian(effective_component(bloch, n), effective_component(symbolic, n))
  end
  @test same_liouvillian(effective_generator(bloch), effective_generator(symbolic))
  @test symbolically_same(newborn_channel(bloch).rate, newborn.rate)
end

@testset "symbolic GKSLNormalForm with some parameters exact" begin
  h = FockSpace(:c)
  a = Destroy(h, :a)
  H = 1 // 2 * a' * a + 3 // 10 * a' * a * a' * a + F * (a + a') * cos(ω * t)
  symbolic = floquet_expansion(H, ω, t, GKSLNormalForm(), 5; channels=(jump(a, κ),))
  newborn = newborn_channel(symbolic)
  @test symbolically_same(newborn.rate, convert(FE.SQA.CNum, 9 // 50 * κ * F^2 / ω^4))
  @test iszero(FE.SQA.simplify(newborn.operator - (a * a - 2 * a' * a)))

  Hn = 1 // 2 * a' * a + 3 // 10 * a' * a * a' * a + 2 // 5 * (a + a') * cos(ω * t)
  exact = floquet_expansion(Hn, ω, t, GKSLNormalForm(), 5; channels=(jump(a, 4 // 5),))
  @test nearly_same(effective_generator(symbolic), effective_generator(exact))
end

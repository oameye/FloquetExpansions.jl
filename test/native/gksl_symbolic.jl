using Test
using Symbolics: @variables
using FloquetExpansions
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions
const SQA = FE.SQA

@variables ω::Real t::Real Δ::Real F::Real κ::Real χ::Real

const QUBIT = PauliSpace(:gksl_symbolic_qubit)
const σx = Pauli(QUBIT, :σ, 1)
const σy = Pauli(QUBIT, :σ, 2)
const σz = Pauli(QUBIT, :σ, 3)
const CAVITY = FockSpace(:gksl_symbolic_cavity)
const a = Destroy(CAVITY, :a)

qubit_hamiltonian(Δ, F) = Δ * σz / 2 + F * cos(ω * t) * σx
decay(κ) = jump((σx - im * σy) / 2, κ)
kerr_hamiltonian(χ, F) = (1 // 2) * a' * a + χ * a' * a * a' * a + F * (a + a') * cos(ω * t)

function evaluated(L::Liouvillian, rules)
  return Dict((string(l), string(r)) => SQA.substitute(c, rules) for (l, r, c) in terms(L))
end

function agree(A::Liouvillian, B::Liouvillian, rules_a, rules_b)
  x, y = evaluated(A, rules_a), evaluated(B, rules_b)
  return keys(x) == keys(y) && all(iszero(SQA.simplify(x[k] - y[k])) for k in keys(x))
end

function agree(A::FloquetExpansion, B::FloquetExpansion, rules_a, rules_b, order::Int)
  effective = all(
    agree(effective_component(A, n), effective_component(B, n), rules_a, rules_b) for
    n in 0:(order - 1)
  )
  micro = all(1:(order - 1)) do n
    x, y = micromotion(A, n).components, micromotion(B, n).components
    return all(
      agree(get(x, k, zero(Liouvillian)), get(y, k, zero(Liouvillian)), rules_a, rules_b)
      for k in union(keys(x), keys(y))
    )
  end
  jumps = length(channels(A)) == length(channels(B))
  rates =
    jumps && all(
      iszero(
        SQA.simplify(SQA.substitute(c.rate, rules_a) - SQA.substitute(d.rate, rules_b))
      ) && iszero(
        SQA.simplify(
          SQA.substitute(c.operator, rules_a) - SQA.substitute(d.operator, rules_b)
        ),
      ) for (c, d) in zip(channels(A), channels(B))
    )
  return effective && micro && rates
end

same(x, y) = iszero(SQA.simplify(convert(SQA.CNum, x) - convert(SQA.CNum, y)))

@testset "symbolic GKSLNormalForm equals the exact expansion at its verification points" begin
  for algorithm in (HoriDeprit(), BlochFeshbach())
    gauge = GKSLNormalForm(; algorithm, verification=4)
    x = floquet_expansion(qubit_hamiltonian(Δ, F), ω, t, gauge, 5; channels=(decay(κ),))
    certificate = factorization(x)
    @test length(certificate.points) == 4
    @test all(==(1), certificate.signs)
    values = [p[findfirst(isequal(F), certificate.parameters)] for p in certificate.points]
    @test maximum(values) / minimum(values) > 50
    for point in certificate.points
      rules = Dict(zip(certificate.parameters, point))
      exact = floquet_expansion(
        qubit_hamiltonian(rules[Δ], rules[F]),
        ω,
        t,
        GKSLNormalForm(; algorithm),
        5;
        channels=(decay(rules[κ]),),
      )
      @test agree(x, exact, merge(rules, Dict(ω => 7)), Dict(ω => 7), 5)
    end
  end
end

@testset "symbolic GKSLNormalForm equals Van Vleck and Gram wherever no static gauge is needed" begin
  rules = Dict(Δ => 1 // 2, F => 2 // 5, κ => 3 // 10, χ => 3 // 10, ω => 7)
  for (H, channel) in
      ((qubit_hamiltonian(Δ, F), decay(κ)), (kerr_hamiltonian(χ, F), jump(a, κ)))
    symbolic = floquet_expansion(
      H, ω, t, GKSLNormalForm(; verification=2), 4; channels=(channel,)
    )
    vanvleck = floquet_expansion(H, ω, t, VanVleck(), 4; channels=(channel,))
    gram = positive_completion(vanvleck, Gram())
    for n in 0:3
      @test agree(
        effective_component(symbolic, n), effective_component(vanvleck, n), rules, rules
      )
    end
    @test agree(effective_generator(symbolic), effective_generator(gram), rules, rules)
  end
end

@testset "symbolic GKSLNormalForm reproduces the order-four channel of driven Kerr" begin
  x = floquet_expansion(
    kerr_hamiltonian(χ, F),
    ω,
    t,
    GKSLNormalForm(; verification=3),
    5;
    channels=(jump(a, 4 // 5),),
  )
  @test length(channels(x)) == 2
  newborn = channels(x)[2]
  rate = 8 // 5 * χ^2 * F^2 / ω^4
  @test same(newborn.rate, rate)
  @test iszero(SQA.simplify(newborn.operator - (a * a - 2 * a' * a)))
  conditions = positivity_conditions(x)
  for target in (χ, F, ω, rate)
    @test any(c -> same(c, target), conditions)
  end
  @test_throws ArgumentError positive_completion(x, Gram())
end

@testset "symbolic GKSLNormalForm reports and samples the declared sign region" begin
  x = floquet_expansion(
    qubit_hamiltonian(Δ, F),
    ω,
    t,
    GKSLNormalForm(; negative=[F], verification=3),
    4;
    channels=(decay(κ),),
  )
  certificate = factorization(x)
  signs = Dict(zip(certificate.parameters, certificate.signs))
  @test signs[F] == -1 && signs[Δ] == 1 && signs[κ] == 1
  @test all(p -> p[findfirst(isequal(F), certificate.parameters)] < 0, certificate.points)
  conditions = positivity_conditions(x)
  @test any(c -> same(c, -F), conditions)
  @test !any(c -> same(c, F), conditions)
  @test any(c -> same(c, ω), conditions)
  for point in certificate.points
    rules = Dict(zip(certificate.parameters, point))
    exact = floquet_expansion(
      qubit_hamiltonian(rules[Δ], rules[F]),
      ω,
      t,
      GKSLNormalForm(),
      4;
      channels=(decay(rules[κ]),),
    )
    @test agree(x, exact, merge(rules, Dict(ω => 7)), Dict(ω => 7), 4)
  end
end

@testset "symbolic GKSLNormalForm rejects what it cannot certify" begin
  channels_ = (decay(κ),)
  gauge(; kwargs...) = GKSLNormalForm(; kwargs...)
  H = qubit_hamiltonian(Δ, F)
  @test_throws ArgumentError floquet_expansion(
    H, ω, t, gauge(; negative=[χ]), 3; channels=channels_
  )
  @test_throws ArgumentError floquet_expansion(
    H, ω, t, gauge(; negative=[ω]), 3; channels=channels_
  )
  @test_throws ArgumentError floquet_expansion(
    H, ω, t, gauge(; negative=[κ], verification=2), 3; channels=channels_
  )
  @test_throws ArgumentError floquet_expansion(
    sqrt(Δ) * σz / 2 + F * cos(ω * t) * σx, ω, t, gauge(), 3; channels=channels_
  )
  @test_throws ArgumentError gauge(; verification=0)
  many = @variables p1 p2 p3 p4 p5 p6 p7
  @test_throws ArgumentError floquet_expansion(
    sum(many) * σz / 2 + F * cos(ω * t) * σx, ω, t, gauge(), 3; channels=channels_
  )
end

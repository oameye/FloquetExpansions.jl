using Test
using Symbolics: Symbolics, @variables
using FloquetExpansions
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions
const SQA = FE.SQA

@variables ω::Real t::Real Δ::Real χ::Real F::Real κ::Real γ::Real Ω::Real α::Complex

const CAVITY = FockSpace(:gksl_symbolic_cavity)
const A = Destroy(CAVITY, :a)

kerr(Δ, χ, F) = Δ * A' * A + χ * A' * A * A' * A + F * (A + A') * cos(ω * t)

function evaluated(L::Liouvillian, rules)
  return Dict((string(l), string(r)) => SQA.substitute(c, rules) for (l, r, c) in terms(L))
end

function agree(symbolic::Liouvillian, numeric::Liouvillian, rules)
  x = evaluated(symbolic, rules)
  y = evaluated(numeric, Dict(ω => rules[ω]))
  return keys(x) == keys(y) && all(iszero(SQA.simplify(x[k] - y[k])) for k in keys(x))
end

function agree(symbolic::FloquetExpansion, numeric::FloquetExpansion, rules)
  N = order(symbolic)
  frequency = Dict(ω => rules[ω])
  effective = all(
    agree(effective_component(symbolic, n), effective_component(numeric, n), rules) for
    n in 0:(N - 1)
  )
  micro = all(1:(N - 1)) do n
    x, y = micromotion(symbolic, n).components, micromotion(numeric, n).components
    return all(
      agree(get(x, k, zero(Liouvillian)), get(y, k, zero(Liouvillian)), rules) for
      k in union(keys(x), keys(y))
    )
  end
  jumps =
    length(channels(symbolic)) == length(channels(numeric)) &&
    all(zip(channels(symbolic), channels(numeric))) do (c, d)
      return iszero(
        SQA.simplify(SQA.substitute(c.rate, rules) - SQA.substitute(d.rate, frequency))
      ) && iszero(
        SQA.simplify(
          SQA.substitute(c.operator, rules) - SQA.substitute(d.operator, frequency)
        ),
      )
    end
  return effective && micro && jumps
end

@testset "symbolic GKSLNormalForm equals the numeric expansion at substituted points" begin
  for (negative, point) in (
    ((), Dict(Δ => 3 // 7, χ => 2 // 9, F => 5 // 11, κ => 7 // 5, ω => 13)),
    ((Δ,), Dict(Δ => -5 // 7, χ => 3 // 10, F => 2 // 5, κ => 4 // 5, ω => 13)),
  )
    symbolic = floquet_expansion(
      kerr(Δ, χ, F), ω, t, GKSLNormalForm(; negative), 5; channels=(jump(A, κ),)
    )
    numeric = floquet_expansion(
      kerr(point[Δ], point[χ], point[F]),
      ω,
      t,
      GKSLNormalForm(),
      5;
      channels=(jump(A, point[κ]),),
    )
    @test agree(symbolic, numeric, point)
  end
  value = 1 // 3 + (2 // 7) * im
  displaced(α, Δ, χ, F) = α * A' + conj(α) * A + kerr(Δ, χ, F)
  re, ie = (Symbolics.wrap(f(Symbolics.unwrap(α))) for f in (real, imag))
  point = Dict(
    re => 1 // 3, ie => 2 // 7, Δ => 3 // 7, χ => 2 // 9, F => 5 // 11, κ => 7 // 5, ω => 13
  )
  symbolic = floquet_expansion(
    displaced(α, Δ, χ, F), ω, t, GKSLNormalForm(), 4; channels=(jump(A, κ),)
  )
  numeric = floquet_expansion(
    displaced(value, 3 // 7, 2 // 9, 5 // 11),
    ω,
    t,
    GKSLNormalForm(),
    4;
    channels=(jump(A, 7 // 5),),
  )
  @test agree(symbolic, numeric, point)
  qubit = PauliSpace(:gksl_symbolic_qubit)
  σx = Pauli(qubit, :σ, 1)
  σy = Pauli(qubit, :σ, 2)
  σz = Pauli(qubit, :σ, 3)
  σminus = (1 // 2) * (σx - im * σy)
  qubit_model(Δ, F, Ω) = Δ * σz + F * σx + Ω * cos(ω * t) * σz
  point = Dict(Δ => 1 // 2, F => 1 // 5, Ω => 2 // 5, γ => 3 // 10, ω => 13)
  symbolic = floquet_expansion(
    qubit_model(Δ, F, Ω), ω, t, GKSLNormalForm(), 4; channels=(jump(σminus, γ),)
  )
  numeric = floquet_expansion(
    qubit_model(point[Δ], point[F], point[Ω]),
    ω,
    t,
    GKSLNormalForm(),
    4;
    channels=(jump(σminus, point[γ]),),
  )
  @test agree(symbolic, numeric, point)
end

function same(conditions, expected)
  return length(conditions) == length(expected) && all(
    iszero(SQA.simplify(c - convert(SQA.CNum, e))) for (c, e) in zip(conditions, expected)
  )
end

@testset "symbolic GKSLNormalForm states where it holds" begin
  gauge = GKSLNormalForm(; negative=(Δ,))
  @test gauge == GKSLNormalForm(; negative=(Δ,))
  @test hash(FE.ParameterScalar(3 // 7)) == hash(3 // 7)
  loss = floquet_expansion(kerr(Δ, χ, F), ω, t, gauge, 5; channels=(jump(A, κ),))
  @test same(positivity_conditions(loss), (F, -Δ, κ, χ, ω))
  @test isempty(regularity_conditions(loss))
  dephased = floquet_expansion(
    kerr(Δ, χ, F), ω, t, GKSLNormalForm(), 4; channels=(jump(A, κ), jump(A' * A, γ))
  )
  @test same(regularity_conditions(dephased), (κ,))
end

@testset "symbolic GKSLNormalForm rejects what it cannot decide" begin
  H = kerr(Δ, χ, F)
  @test_throws "the sign of γ - κ is not fixed" floquet_expansion(
    H, ω, t, GKSLNormalForm(), 5; channels=(jump(A, κ - γ),)
  )
  @test_throws "does not occur in the generator" floquet_expansion(
    H, ω, t, GKSLNormalForm(; negative=(γ,)), 3; channels=(jump(A, κ),)
  )
  @test_throws "drive frequency is positive" floquet_expansion(
    H, ω, t, GKSLNormalForm(; negative=(ω,)), 3; channels=(jump(A, κ),)
  )
  @test_throws "drive frequency ω appears in a coefficient" floquet_expansion(
    (1 // 2) * ω * A' * A + F * (A + A') * cos(ω * t),
    ω,
    t,
    GKSLNormalForm(),
    3;
    channels=(jump(A, κ),),
  )
  @test_throws "integer powers of a parameter" floquet_expansion(
    H, ω, t, GKSLNormalForm(), 3; channels=(jump(A, sqrt(κ)),)
  )
  @variables p[1:9]::Real
  crowded = sum(p[i] * (A' * A)^i for i in 1:9) + F * (A + A') * cos(ω * t)
  @test_throws "at most 8 symbolic parameters" floquet_expansion(
    crowded, ω, t, GKSLNormalForm(), 3; channels=(jump(A, κ),)
  )
end

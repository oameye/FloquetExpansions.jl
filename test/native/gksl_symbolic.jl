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

lindbladian(H, rate) = hamiltonian_action(H) + rate * dissipator(A)

channel_map(entry) = SQA.SuperOperator(last(entry) * dissipator(SQA.jump(first(entry))))

function evaluated(L::SQA.SuperOperator, rules)
  return Dict(
    (string(SQA.left_operator(term)), string(SQA.right_operator(term))) =>
      SQA.substitute(c, rules) for (term, c) in L
  )
end

function agree(symbolic::SQA.SuperOperator, numeric::SQA.SuperOperator, rules)
  x = evaluated(symbolic, rules)
  y = evaluated(numeric, Dict(ω => rules[ω]))
  vanishing = convert(SQA.CNum, 0)
  return all(
    iszero(SQA.simplify(get(x, k, vanishing) - get(y, k, vanishing))) for
    k in union(keys(x), keys(y))
  )
end

function agree(symbolic::FloquetExpansion, numeric::FloquetExpansion, rules)
  N = order(symbolic)
  effective = all(
    agree(effective_component(symbolic, n), effective_component(numeric, n), rules) for
    n in 0:(N - 1)
  )
  micro = all(1:(N - 1)) do n
    x, y = micromotion(symbolic, n).components, micromotion(numeric, n).components
    zero_map = zero(SQA.SuperOperator)
    return all(
      agree(get(x, k, zero_map), get(y, k, zero_map), rules) for
      k in union(keys(x), keys(y))
    )
  end
  jumps =
    length(channels(symbolic)) == length(channels(numeric)) && all(channels(symbolic)) do c
      matches = d -> agree(channel_map(c), channel_map(d), rules)
      return count(matches, channels(numeric)) == 1
    end
  return effective && micro && jumps
end

@testset "symbolic GKSLNormalForm equals the numeric expansion at substituted points" begin
  for (negative, point) in (
    ((), Dict(Δ => 3 // 7, χ => 2 // 9, F => 5 // 11, κ => 7 // 5, ω => 13)),
    ((Δ,), Dict(Δ => -5 // 7, χ => 3 // 10, F => 2 // 5, κ => 4 // 5, ω => 13)),
  )
    symbolic = floquet_expansion(
      lindbladian(kerr(Δ, χ, F), κ), ω, t, GKSLNormalForm(; negative), 5
    )
    numeric = floquet_expansion(
      lindbladian(kerr(point[Δ], point[χ], point[F]), point[κ]), ω, t, GKSLNormalForm(), 5
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
    lindbladian(displaced(α, Δ, χ, F), κ), ω, t, GKSLNormalForm(), 4
  )
  numeric = floquet_expansion(
    lindbladian(displaced(value, 3 // 7, 2 // 9, 5 // 11), 7 // 5),
    ω,
    t,
    GKSLNormalForm(),
    4,
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
    hamiltonian_action(qubit_model(Δ, F, Ω)) + γ * dissipator(σminus),
    ω,
    t,
    GKSLNormalForm(),
    4,
  )
  numeric = floquet_expansion(
    hamiltonian_action(qubit_model(point[Δ], point[F], point[Ω])) +
    point[γ] * dissipator(σminus),
    ω,
    t,
    GKSLNormalForm(),
    4,
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
  loss = floquet_expansion(lindbladian(kerr(Δ, χ, F), κ), ω, t, gauge, 5)
  @test same(positivity_conditions(loss), (F, -Δ, κ, χ, ω))
  @test isempty(regularity_conditions(loss))
  dephased = floquet_expansion(
    lindbladian(kerr(Δ, χ, F), κ) + γ * dissipator(A' * A), ω, t, GKSLNormalForm(), 4
  )
  @test same(regularity_conditions(dephased), (κ,))
end

@testset "symbolic GKSLNormalForm rejects what it cannot decide" begin
  H = kerr(Δ, χ, F)
  @test_throws "the sign of γ - κ is not fixed" floquet_expansion(
    lindbladian(H, κ - γ), ω, t, GKSLNormalForm(), 5
  )
  @test_throws "does not occur in the generator" floquet_expansion(
    lindbladian(H, κ), ω, t, GKSLNormalForm(; negative=(γ,)), 3
  )
  @test_throws "drive frequency is positive" floquet_expansion(
    lindbladian(H, κ), ω, t, GKSLNormalForm(; negative=(ω,)), 3
  )
  @test_throws "drive frequency ω appears in a coefficient" floquet_expansion(
    lindbladian((1 // 2) * ω * A' * A + F * (A + A') * cos(ω * t), κ),
    ω,
    t,
    GKSLNormalForm(),
    3,
  )
  @test_throws "integer powers of a parameter" floquet_expansion(
    lindbladian(H, sqrt(κ)), ω, t, GKSLNormalForm(), 3
  )
  @variables p[1:9]::Real
  crowded = sum(p[i] * (A' * A)^i for i in 1:9) + F * (A + A') * cos(ω * t)
  @test_throws "at most 8 symbolic parameters" floquet_expansion(
    lindbladian(crowded, κ), ω, t, GKSLNormalForm(), 3
  )
end

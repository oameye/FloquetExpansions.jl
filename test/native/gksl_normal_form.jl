using Test
using LinearAlgebra: Hermitian, eigmin
using Symbolics: @variables, substitute
using FloquetExpansions
using FloquetExpansions: FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra

const FE = FloquetExpansions
const SQA = SecondQuantizedAlgebra

@variables ω::Real t::Real

include(joinpath(@__DIR__, "..", "helpers", "shared.jl"))

function has_channel(jumps, expected::Lindbladian)
  return count(entry -> same(channel_map(entry), SQA.SuperOperator(expected)), jumps) == 1
end

function kerr_hamiltonian()
  h = FockSpace(:c)
  a = Destroy(h, :a)
  H = (1 // 2) * a' * a + (3 // 10) * a' * a * a' * a + (2 // 5) * (a + a') * cos(ω * t)
  return H, a
end

function graded_product(A, B)
  C = [Dict{Int,SQA.SuperOperator}() for _ in A]
  for i in eachindex(A), j in 1:(length(A) - i), (k, a) in A[i], (l, b) in B[j]
    C[i + j][k + l] = get(C[i + j], k + l, zero(a)) + SQA.compose(a, b)
  end
  return C
end

graded_sum(A, B) = [mergewith(+, a, b) for (a, b) in zip(A, B)]

function graded_exp_minus_one(K)
  result = K
  power = K
  for p in 2:length(K)
    power = graded_product(power, K)
    scaled = [Dict(k => (1 // factorial(p)) * v for (k, v) in d) for d in power]
    result = graded_sum(result, scaled)
  end
  return result
end

function micromotion_series(expansion, sign)
  return [
    Dict(k => sign * v for (k, v) in micromotion(expansion, n).components) for
    n in 1:(expansion.order - 1)
  ]
end

@testset "GKSLNormalForm is a static similarity of VanVleck where the static gauge is nonzero" begin
  qubit = NLevelSpace(:static_gauge_qubit, 2)
  σ(i, j) = Transition(qubit, :σ, i, j)
  σz = σ(1, 1) - σ(2, 2)
  σx = σ(1, 2) + σ(2, 1)
  σy = -im * σ(1, 2) + im * σ(2, 1)
  H = (1 // 2) * σz + cos(ω * t) * σx + (1 // 3) * cos(2 * ω * t) * σy
  L = hamiltonian_action(H) + (1 // 5) * dissipator(σ(1, 2)) + (1 // 7) * dissipator(σz)
  N = 4
  for algorithm in (HoriDeprit(), BlochFeshbach())
    native = floquet_expansion(L, ω, t, GKSLNormalForm(; algorithm), N)
    reference = floquet_expansion(L, ω, t, VanVleck(; algorithm), N)
    A = graded_exp_minus_one(micromotion_series(reference, -1))
    B = graded_exp_minus_one(micromotion_series(native, 1))
    W = graded_sum(graded_sum(A, B), graded_product(A, B))
    for n in 1:(N - 1), (k, value) in W[n]
      k == 0 || @test same(value, zero(value))
    end
    S = [get(W[n], 0, zero(SQA.SuperOperator)) for n in 1:(N - 1)]
    @test !same(S[2], zero(SQA.SuperOperator))
    E(x, n) = effective_component(x, n)
    for n in 1:(N - 1)
      @test same(
        E(native, n) + sum(SQA.compose(S[j], E(native, n - j)) for j in 1:n),
        E(reference, n) + sum(SQA.compose(E(reference, n - j), S[j]) for j in 1:n),
      )
    end
  end
end

@testset "GKSLNormalForm coincides with VanVleck wherever no static gauge is needed" begin
  H, a = kerr_hamiltonian()
  for algorithm in (HoriDeprit(), BlochFeshbach()),
    generator in (H, hamiltonian_action(H) + (4 // 5) * dissipator(a))

    native = floquet_expansion(generator, ω, t, GKSLNormalForm(; algorithm), 4)
    reference = floquet_expansion(generator, ω, t, VanVleck(; algorithm), 4)
    for n in 0:3
      @test same(effective_component(native, n), effective_component(reference, n))
    end
    for n in 1:3
      m, r = micromotion(native, n), micromotion(reference, n)
      for k in union(keys(m.components), keys(r.components))
        @test same(m[k], r[k])
      end
    end
  end
end

@testset "GKSLNormalForm driven Kerr is GKSL with a drive-induced channel" begin
  H, a = kerr_hamiltonian()
  L = hamiltonian_action(H) + (4 // 5) * dissipator(a)
  native = floquet_expansion(L, ω, t, GKSLNormalForm(), 5)
  jumps = channels(native)
  @test length(jumps) == 2
  @test has_channel(jumps, (72 // 3125) / ω^4 * dissipator(a * a - 2 * a' * a))
  @test same(effective_generator(native), SQA.SuperOperator(Lindbladian(native)))
  C = kossakowski(native)
  for value in (3.0, 10.0, 100.0)
    numeric = [ComplexF64(SQA.to_complex(SQA.substitute(c, Dict(ω => value)))) for c in C]
    @test eigmin(Hermitian(numeric)) >= -1e-12
  end
  @test_throws ArgumentError positive_completion(native, Gram())
  bf = floquet_expansion(L, ω, t, GKSLNormalForm(; algorithm=BlochFeshbach()), 5)
  @test all(same(effective_component(bf, n), effective_component(native, n)) for n in 0:4)
end

@testset "GKSLNormalForm generates one-photon loss from driven two-photon loss" begin
  @variables Δ::Real K::Real ε::Real κ::Real
  h = FockSpace(:two_photon)
  a = Destroy(h, :a)
  H = Δ * a' * a + K * a' * a' * a * a + ε * cos(ω * t) * (a + a')
  L = hamiltonian_action(H) + κ * dissipator(a * a)
  native = floquet_expansion(L, ω, t, GKSLNormalForm(), 3)
  jumps = channels(native)
  @test length(jumps) == 2
  @test has_channel(jumps, κ * dissipator(a * a))
  @test has_channel(jumps, (2 * ε^2 * κ / ω^2) * dissipator(a))
end

@testset "GKSLNormalForm on a driven oscillator in position and momentum" begin
  phase = PhaseSpace(:q)
  x = Position(phase, :x)
  p = Momentum(phase, :p)
  H0 = (1 // 2) * p * p + (1 // 2) * x * x
  H = H0 + (1 // 5) * x * cos(ω * t)
  loss = (1 // 5) * dissipator(x + im * p)
  native = floquet_expansion(hamiltonian_action(H) + loss, ω, t, GKSLNormalForm(), 3)
  reference = floquet_expansion(hamiltonian_action(H) + loss, ω, t, VanVleck(), 3)
  @test same(
    effective_component(native, 0), SQA.SuperOperator(hamiltonian_action(H0) + loss)
  )
  @test all(
    same(effective_component(native, n), effective_component(reference, n)) for n in 0:2
  )
  @test same(effective_generator(native), SQA.SuperOperator(Lindbladian(native)))
  C = kossakowski(native)
  for value in (3.0, 10.0)
    numeric = [ComplexF64(SQA.to_complex(SQA.substitute(c, Dict(ω => value)))) for c in C]
    @test eigmin(Hermitian(numeric)) >= -1e-12
  end
end

@testset "GKSLNormalForm rejects what it cannot expand exactly" begin
  H, a = kerr_hamiltonian()
  loss = (4 // 5) * dissipator(a)
  @test_throws "order must be >= 1" floquet_expansion(
    hamiltonian_action(H) + loss, ω, t, GKSLNormalForm(), 0
  )
  @test_throws "no static gauge of polynomial degree" floquet_expansion(
    hamiltonian_action(H) + loss, ω, t, GKSLNormalForm(), 7
  )
  @test_throws "floating-point" floquet_expansion(
    hamiltonian_action(H + 0.4 * (a + a') * cos(ω * t)) + loss, ω, t, GKSLNormalForm(), 3
  )
  @test_throws "floating-point" floquet_expansion(
    hamiltonian_action(H) + 0.8 * dissipator(a), ω, t, GKSLNormalForm(), 3
  )
  shared = FockSpace(:shared)
  a = Destroy(shared, :a)
  b = Destroy(shared, :b)
  two_names =
    (1 // 2) * a' * a + (3 // 10) * b' * b * b' * b + (2 // 5) * (a + a') * cos(ω * t)
  @test_throws "single operator name" floquet_expansion(
    hamiltonian_action(two_names) + loss, ω, t, GKSLNormalForm(), 3
  )
  not_hermitian = PeriodicGenerator(
    Dict(
      0 => SQA.SuperOperator(hamiltonian_action(a)),
      1 => SQA.SuperOperator(hamiltonian_action(a + a')),
      -1 => SQA.SuperOperator(hamiltonian_action(a + a')),
    ),
    ω,
  )
  @test_throws "not of GKSL form" floquet_expansion(not_hermitian, GKSLNormalForm(), 3)
  hc = FockSpace(:cavity) ⊗ NLevelSpace(:atoms, 2, 1)
  @variables N::Int
  i = Index(hc, :i, N, NLevelSpace(:atoms, 2, 1))
  c = Destroy(hc, :c, 1)
  σ(j, k) = IndexedOperator(Transition(hc, :σ, j, k, 2), i)
  collective = Σ(cos(ω * t) * (c' * σ(1, 2) + c * σ(2, 1)), i)
  ensemble = SQA.liouvillian(c' * c + collective, [c]; rates=[1 // 5])
  @test_throws "symbolic index sums" floquet_expansion(ensemble, ω, t, GKSLNormalForm(), 2)
  @test_throws "symbolic index sums" floquet_expansion(
    period_map(ensemble, ω, t), VanVleck(), 1
  )
  single = SQA.liouvillian(cos(ω * t) * (c' * σ(1, 2) + c * σ(2, 1)), [c]; rates=[1 // 5])
  @test_throws "indexed operators" floquet_expansion(single, ω, t, GKSLNormalForm(), 2)
end

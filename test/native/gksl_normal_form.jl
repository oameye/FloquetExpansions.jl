using Test
using LinearAlgebra: Hermitian, eigmin
using Symbolics: @variables, substitute
using FloquetExpansions
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions

@variables ω::Real t::Real

function same(A::Liouvillian, B::Liouvillian)
  return FE.liouvillian_iszero(FE.canonical_liouvillian(FE.SQA.simplify(A - B)))
end
same(A, B) = iszero(FE.SQA.simplify(FE.qadd(A) - FE.qadd(B)))

function kerr_hamiltonian()
  h = FockSpace(:c)
  a = Destroy(h, :a)
  H = (1 // 2) * a' * a + (3 // 10) * a' * a * a' * a + (2 // 5) * (a + a') * cos(ω * t)
  return H, a
end

function graded_product(A, B)
  C = [Dict{Int,Liouvillian}() for _ in A]
  for i in eachindex(A), j in 1:(length(A) - i), (k, a) in A[i], (l, b) in B[j]
    C[i + j][k + l] = get(C[i + j], k + l, zero(a)) + compose(a, b)
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
  jumps = (jump(σ(1, 2), 1 // 5), jump(σz, 1 // 7))
  N = 4
  for algorithm in (HoriDeprit(), BlochFeshbach())
    native = floquet_expansion(H, ω, t, GKSLNormalForm(; algorithm), N; channels=jumps)
    reference = floquet_expansion(H, ω, t, VanVleck(; algorithm), N; channels=jumps)
    A = graded_exp_minus_one(micromotion_series(reference, -1))
    B = graded_exp_minus_one(micromotion_series(native, 1))
    W = graded_sum(graded_sum(A, B), graded_product(A, B))
    for n in 1:(N - 1), (k, value) in W[n]
      k == 0 || @test same(value, zero(value))
    end
    S = [get(W[n], 0, zero(Liouvillian)) for n in 1:(N - 1)]
    @test !same(S[2], zero(Liouvillian))
    E(x, n) = effective_component(x, n)
    for n in 1:(N - 1)
      @test same(
        E(native, n) + sum(compose(S[j], E(native, n - j)) for j in 1:n),
        E(reference, n) + sum(compose(E(reference, n - j), S[j]) for j in 1:n),
      )
    end
  end
end

@testset "GKSLNormalForm coincides with VanVleck wherever no static gauge is needed" begin
  H, a = kerr_hamiltonian()
  for algorithm in (HoriDeprit(), BlochFeshbach()), channels in ((), (jump(a, 4 // 5),))
    native = floquet_expansion(H, ω, t, GKSLNormalForm(; algorithm), 4; channels)
    reference = floquet_expansion(H, ω, t, VanVleck(; algorithm), 4; channels)
    for n in 0:3
      @test same(effective_component(native, n), effective_component(reference, n))
    end
    for n in 1:3
      m, r = micromotion(native, n), micromotion(reference, n)
      for k in union(keys(m.components), keys(r.components))
        @test same(
          get(m.components, k, zero(m.zero_component)),
          get(r.components, k, zero(r.zero_component)),
        )
      end
    end
  end
end

@testset "GKSLNormalForm driven Kerr is GKSL with a drive-induced channel" begin
  H, a = kerr_hamiltonian()
  native = floquet_expansion(H, ω, t, GKSLNormalForm(), 5; channels=(jump(a, 4 // 5),))
  jumps = channels(native)
  @test length(jumps) == 2
  @test same(jumps[2].operator, a * a - 2 * a' * a)
  @test isequal(FE.SQA.simplify(jumps[2].rate * ω^4), 72 // 3125)
  @test same(
    effective_generator(native), liouvillian(hamiltonian(native); channels=Tuple(jumps))
  )
  C = kossakowski(native)
  for value in (3.0, 10.0, 100.0)
    numeric = [
      ComplexF64(FE.SQA.to_complex(FE.SQA.substitute(c, Dict(ω => value)))) for c in C
    ]
    @test eigmin(Hermitian(numeric)) >= -1e-12
  end
  @test_throws ArgumentError positive_completion(native, Gram())
  bf = floquet_expansion(
    H, ω, t, GKSLNormalForm(; algorithm=BlochFeshbach()), 5; channels=(jump(a, 4 // 5),)
  )
  @test all(same(effective_component(bf, n), effective_component(native, n)) for n in 0:4)
end

@testset "GKSLNormalForm on a driven oscillator in position and momentum" begin
  phase = PhaseSpace(:q)
  x = Position(phase, :q)
  p = Momentum(phase, :q)
  H0 = (1 // 2) * p * p + (1 // 2) * x * x
  H = H0 + (1 // 5) * x * cos(ω * t)
  jumps = (jump(x + im * p, 1 // 5),)
  native = floquet_expansion(H, ω, t, GKSLNormalForm(), 3; channels=jumps)
  reference = floquet_expansion(H, ω, t, VanVleck(), 3; channels=jumps)
  @test same(effective_component(native, 0), liouvillian(H0; channels=jumps))
  @test all(
    same(effective_component(native, n), effective_component(reference, n)) for n in 0:2
  )
  @test same(
    effective_generator(native),
    liouvillian(hamiltonian(native); channels=Tuple(channels(native))),
  )
  C = kossakowski(native)
  for value in (3.0, 10.0)
    numeric = [
      ComplexF64(FE.SQA.to_complex(FE.SQA.substitute(c, Dict(ω => value)))) for c in C
    ]
    @test eigmin(Hermitian(numeric)) >= -1e-12
  end
end

@testset "GKSLNormalForm rejects what it cannot expand exactly" begin
  H, a = kerr_hamiltonian()
  jumps = (jump(a, 4 // 5),)
  @test_throws "order must be >= 1" floquet_expansion(
    H, ω, t, GKSLNormalForm(), 0; channels=jumps
  )
  @test_throws "no static gauge of polynomial degree" floquet_expansion(
    H, ω, t, GKSLNormalForm(), 7; channels=jumps
  )
  @test_throws "floating-point" floquet_expansion(
    H + 0.4 * (a + a') * cos(ω * t), ω, t, GKSLNormalForm(), 3; channels=jumps
  )
  @test_throws "floating-point" floquet_expansion(
    H, ω, t, GKSLNormalForm(), 3; channels=(jump(a, 0.8),)
  )
  shared = FockSpace(:shared)
  a = Destroy(shared, :a)
  b = Destroy(shared, :b)
  two_names =
    (1 // 2) * a' * a + (3 // 10) * b' * b * b' * b + (2 // 5) * (a + a') * cos(ω * t)
  @test_throws "single operator name" floquet_expansion(
    two_names, ω, t, GKSLNormalForm(), 3; channels=jumps
  )
  not_hermitian = PeriodicGenerator(
    Dict(
      0 => hamiltonian_action(a),
      1 => hamiltonian_action(a + a'),
      -1 => hamiltonian_action(a + a'),
    ),
    ω,
  )
  @test_throws "not of GKSL form" floquet_expansion(not_hermitian, GKSLNormalForm(), 3)
end

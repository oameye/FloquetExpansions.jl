using Test
using FloquetExpansions
using Symbolics: Symbolics, @variables
using LinearAlgebra: I, kron

include(joinpath(@__DIR__, "helpers", "shared.jl"))

space = FockSpace(:liouvillian)
a = Destroy(space, :a)
@variables γ::Real ω::Real t::Real

function superoperator_matrix(L::Liouvillian, space, d::Int, substitutions::AbstractDict)
  matrix = zeros(ComplexF64, d^2, d^2)
  for j in 1:d, i in 1:d
    rho = Transition(space, :sigma, i, j)
    matrix[:, (j - 1) * d + i] = vec(tomatrix(L(rho), d, substitutions))
  end
  return matrix
end

function independent_dissipator_matrix(C::Matrix{ComplexF64})
  norm = C' * C
  identity = Matrix{ComplexF64}(I, size(C, 1), size(C, 2))
  return kron(conj(C), C) - 0.5 * (kron(identity, norm) + kron(transpose(norm), identity))
end

function independent_liouvillian_matrix(H, C, J, rate)
  identity = Matrix{ComplexF64}(I, size(H, 1), size(H, 2))
  coherent = -im * (kron(identity, H) - kron(transpose(H), identity))
  return coherent +
         independent_dissipator_matrix(C) +
         rate * independent_dissipator_matrix(J)
end

@testset "Liouvillian channel constructors" begin
  H = a' * a
  coherent = hamiltonian_action(H)
  complete = dissipator(a)
  weighted = γ * dissipator(a')

  L = Liouvillian(H; channels=(collapse(a), jump(a'; rate=γ)))

  @test L == coherent + complete + weighted
  @test L isa Liouvillian
  @test @inferred(hamiltonian_action(H)) isa Liouvillian
  @test @inferred(dissipator(a)) isa Liouvillian
  @test @inferred(γ * dissipator(a')) isa Liouvillian
  @test @inferred(Liouvillian(H; channels=(jump(a; rate=γ),))) isa Liouvillian
end

@testset "Liouvillian arithmetic collects equal actions" begin
  L = hamiltonian_action(a' * a)

  @test iszero(L - L)
  @test zero(L) + L == L
  @test 2 * L == L + L
  @test L * 2 == L + L
  @test SQA.simplify(L - L) == zero(L)
end

@testset "zero operator factors produce zero maps" begin
  zero_operator = zero(SQA.QAdd)
  @test iszero(hamiltonian_action(zero_operator))
  @test iszero(dissipator(zero_operator))
  @test iszero(compose(hamiltonian_action(a), hamiltonian_action(zero_operator)))
end

@testset "Liouvillian composition is map composition" begin
  coherent = hamiltonian_action(a' * a)
  dissipative = dissipator(a)
  composed = compose(coherent, dissipative)

  @test composed isa Liouvillian
  @test !iszero(composed)
  @test compose(coherent, zero(dissipative)) == zero(coherent)
  @test compose(zero(coherent), dissipative) == zero(coherent)
  @test SQA.commutator(coherent, dissipative) ==
    compose(coherent, dissipative) - compose(dissipative, coherent)
  @test compose(coherent, dissipative) != compose(dissipative, coherent)
end

@testset "Liouvillian channel adapters" begin
  @test Liouvillian(a; channels=(collapse(a),)) == hamiltonian_action(a) + dissipator(a)
  @test Liouvillian(a; channels=(jump(a; rate=γ),)) ==
    hamiltonian_action(a) + γ * dissipator(a)
  @test Liouvillian(zero(SQA.QAdd); channels=(collapse(2a),)) == dissipator(2a)
  @test_throws UndefKeywordError jump(a)
  @test_throws MethodError Liouvillian(a; channels=(a,))
end

@testset "Liouvillians use the common van Vleck expansion" begin
  static = Liouvillian(a' * a; channels=(jump(a; rate=γ),))
  driven = Liouvillian(a)
  generator = PeriodicGenerator(Dict(0 => static, 1 => driven, -1 => driven), ω)
  expansion = floquet_expansion(generator, VanVleck(), 2)

  @test expansion isa FloquetExpansion
  @test effective_generator(expansion) isa Liouvillian
  @test micromotion(expansion) isa PeriodicGenerator{Liouvillian}
  @test effective_generator(expansion, 0) == static
end

@testset "periodic channel forms lower through the public Liouvillian seam" begin
  H = a' * a + cos(ω * t) * (a + a')
  collapse_operator = expim(-ω * t) * a
  jump_operator = expim(-ω * t) * a
  rate = 1 + cos(ω * t)

  native = Liouvillian(H; channels=(collapse(collapse_operator), jump(jump_operator; rate)))
  periodic = harmonics(native, ω, t)

  @test periodic isa PeriodicGenerator{Liouvillian}
  @test @inferred(harmonics(native, ω, t)) isa PeriodicGenerator{Liouvillian}
  @test harmonics(periodic(t), ω, t) == periodic
  @test periodic[0] == Liouvillian(a' * a; channels=(collapse(a), jump(a; rate=1)))
  expected_oscillatory = hamiltonian_action((1 // 2) * (a + a')) + (1 // 2) * dissipator(a)
  @test periodic[1] == expected_oscillatory
  @test periodic[-1] == expected_oscillatory

  direct = floquet_expansion(native, ω, t, VanVleck(), 1)
  keyword = floquet_expansion(
    H,
    ω,
    t,
    VanVleck(),
    1;
    channels=(collapse(collapse_operator), jump(jump_operator; rate)),
  )
  explicit = floquet_expansion(periodic, VanVleck(), 1)

  @test direct isa FloquetExpansion
  @test effective_generator(direct) == effective_generator(explicit)
  @test effective_generator(keyword) == effective_generator(explicit)

  complete_only = harmonics(
    Liouvillian(zero(SQA.QAdd); channels=(collapse(collapse_operator),)), ω, t
  )
  weighted_only = harmonics(
    Liouvillian(zero(SQA.QAdd); channels=(jump(jump_operator; rate=1),)), ω, t
  )
  @test complete_only == weighted_only
end

@testset "operator and rate harmonics convolve in periodic channels" begin
  collapse_operator = a + cos(ω * t) * a'
  jump_operator = a + expim(-ω * t) * a'
  rate = 1 + cos(ω * t)
  periodic = harmonics(
    Liouvillian(
      zero(SQA.QAdd); channels=(collapse(collapse_operator), jump(jump_operator; rate))
    ),
    ω,
    t,
  )

  @test sort!(collect(keys(periodic))) == [-2, -1, 0, 1, 2]
  @test harmonics(periodic(t), ω, t) == periodic
end

@testset "a zero periodic Liouvillian keeps its public zero prototype" begin
  periodic = harmonics(zero(Liouvillian), ω, t)

  @test periodic isa PeriodicGenerator{Liouvillian}
  @test iszero(periodic)
  @test time_average(periodic) == zero(Liouvillian)
  @test iszero(
    floquet_expansion(zero(Liouvillian), ω, t, VanVleck(), 1) |> effective_generator
  )
end

@testset "periodic Liouvillian action agrees with an independent matrix model" begin
  finite_space = NLevelSpace(:finite_validation, 2)
  sigma = Transition(finite_space, :sigma, 1, 2)
  H = sigma + sigma'
  collapse_operator = expim(-ω * t) * sigma
  jump_operator = sigma + expim(ω * t) * sigma'
  rate = 1 + expim(ω * t) + expim(-ω * t)
  native = Liouvillian(H; channels=(collapse(collapse_operator), jump(jump_operator; rate)))
  periodic = harmonics(native, ω, t)
  substitutions = Dict(ω => 3.0, t => 0.37)

  H_matrix = tomatrix(H, 2, substitutions)
  C_matrix = tomatrix(collapse_operator, 2, substitutions)
  J_matrix = tomatrix(jump_operator, 2, substitutions)
  rate_value = real(1 + cis(3.0 * 0.37) + cis(-3.0 * 0.37))
  expected = independent_liouvillian_matrix(H_matrix, C_matrix, J_matrix, rate_value)
  native_matrix = superoperator_matrix(native, finite_space, 2, substitutions)
  periodic_matrix = superoperator_matrix(periodic(t), finite_space, 2, substitutions)

  @test native_matrix ≈ expected
  @test periodic_matrix ≈ expected
end

using Test
using FloquetExpansions
using Symbolics: Num, Symbolics, @variables
using LinearAlgebra: I, kron

include(joinpath(@__DIR__, "..", "helpers", "shared.jl"))

space = FockSpace(:superoperator)
a = Destroy(space, :a)
@variables γ::Real δ::Real ω::Real t::Real t₀::Real

function coeff_to_complex(coeff, substitutions::AbstractDict)
  dummy = coeff * one(SQA.QAdd)
  return tomatrix(dummy, 1, substitutions)[1, 1]
end

function superoperator_matrix(S::SQA.SuperOperator, d::Int, substitutions::AbstractDict)
  matrix = zeros(ComplexF64, d^2, d^2)
  for (term, coeff) in S
    left_mat = tomatrix(SQA.left_operator(term), d, substitutions)
    right_mat = tomatrix(SQA.right_operator(term), d, substitutions)
    c = coeff_to_complex(coeff, substitutions)
    iszero(c) && continue
    matrix .+= c .* kron(transpose(right_mat), left_mat)
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

maps_agree(A, B) = iszero(SQA.expand_completeness(A - B))

@testset "Lindbladian channel rates are real and nonnegative by physical assumption" begin
  expand(L) = floquet_expansion(L, ω, t, VanVleck(), 1)
  phase_rate = 1 + expim(ω * t) + expim(-ω * t)

  for rate in (2, 2 + 0im, γ, 1 + cos(ω * t), phase_rate, -γ)
    @test expand(rate * dissipator(a)) isa FloquetExpansion
  end
  @test sort!(collect(keys(harmonics(phase_rate * dissipator(a), ω, t)))) == [-1, 0, 1]

  for rate in (-1, -1.0, -2 + 0im, Num(-1))
    @test_throws "channel rate must be nonnegative" expand(rate * dissipator(a))
  end
  for rate in (im, γ + im * δ, expim(ω * t₀))
    @test_throws "channel rate must be provably real" expand(rate * dissipator(a))
  end
end

@testset "Lindbladian ingest rejects correlated channels and non-Hermitian drives" begin
  @test_throws ArgumentError floquet_expansion(dissipator(a, a'), ω, t, VanVleck(), 1)
  @test_throws ArgumentError floquet_expansion(
    hamiltonian_action(im * (a' * a)), ω, t, VanVleck(), 1
  )
  @test floquet_expansion(SQA.SuperOperator(dissipator(a, a')), ω, t, VanVleck(), 1) isa
    FloquetExpansion
end

@testset "Kossakowski coordinates read the normalized channel rate" begin
  L = dissipator(2a)
  frame = JumpBasis(a)

  D = @inferred kossakowski(L, frame)

  @test D[1, 1] == 4
  @test iszero(hamiltonian(Lindbladian(SQA.SuperOperator(L), frame)))
end

@testset "superoperator harmonics apply finite-level completeness" begin
  q = NLevelSpace(:atom, 2)
  σ11 = Transition(q, :σ, 1, 1)
  σ22 = Transition(q, :σ, 2, 2)
  σ12 = Transition(q, :σ, 1, 2)
  identity_commutator = SQA.spre(σ11 + σ22) - SQA.spost(σ11 + σ22)

  @test iszero(harmonics(identity_commutator, ω, t))
  @test harmonics(SQA.SuperOperator(dissipator(σ12)) + identity_commutator, ω, t) ==
    harmonics(dissipator(σ12), ω, t)
end

@testset "superoperators use the common van Vleck expansion" begin
  static = SQA.SuperOperator(hamiltonian_action(a' * a) + γ * dissipator(a))
  driven = SQA.SuperOperator(hamiltonian_action(a))
  generator = PeriodicGenerator(Dict(0 => static, 1 => driven, -1 => driven), ω)
  expansion = floquet_expansion(generator, VanVleck(), 2)

  @test expansion isa FloquetExpansion
  @test effective_generator(expansion) isa SQA.SuperOperator
  @test micromotion(expansion) isa PeriodicGenerator{SQA.SuperOperator}
  @test maps_agree(effective_component(expansion, 0), static)
end

@testset "periodic Lindbladians lower through the public superoperator seam" begin
  H = a' * a + cos(ω * t) * (a + a')
  collapse_operator = expim(-ω * t) * a
  jump_operator = expim(-ω * t) * a
  rate = 1 + cos(ω * t)

  native =
    hamiltonian_action(H) + dissipator(collapse_operator) + rate * dissipator(jump_operator)
  periodic = harmonics(native, ω, t)

  @test periodic isa PeriodicGenerator{SQA.SuperOperator}
  @test @inferred(harmonics(native, ω, t)) isa PeriodicGenerator{SQA.SuperOperator}
  @test harmonics(periodic(t), ω, t) == periodic
  @test maps_agree(
    periodic[0], SQA.SuperOperator(hamiltonian_action(a' * a) + 2 * dissipator(a))
  )
  expected_oscillatory = SQA.SuperOperator(
    hamiltonian_action((1 // 2) * (a + a')) + (1 // 2) * dissipator(a)
  )
  @test maps_agree(periodic[1], expected_oscillatory)
  @test maps_agree(periodic[-1], expected_oscillatory)

  direct = floquet_expansion(native, ω, t, VanVleck(), 1)
  lowered = floquet_expansion(SQA.SuperOperator(native), ω, t, VanVleck(), 1)
  explicit = floquet_expansion(periodic, VanVleck(), 1)

  @test direct isa FloquetExpansion
  @test effective_generator(direct) == effective_generator(explicit)
  @test effective_generator(lowered) == effective_generator(explicit)
end

@testset "operator and rate harmonics convolve in periodic channels" begin
  collapse_operator = a + cos(ω * t) * a'
  jump_operator = a + expim(-ω * t) * a'
  rate = 1 + cos(ω * t)
  periodic = harmonics(
    dissipator(collapse_operator) + rate * dissipator(jump_operator), ω, t
  )

  @test sort!(collect(keys(periodic))) == [-2, -1, 0, 1, 2]
  @test harmonics(periodic(t), ω, t) == periodic
end

@testset "general superoperator coefficients keep the shared Fourier convention" begin
  coefficient = expim(ω * t₀) * cos(ω * t)
  periodic = harmonics(coefficient * SQA.SuperOperator(dissipator(a)), ω, t)
  expected = (1 // 2) * expim(ω * t₀) * SQA.SuperOperator(dissipator(a))

  @test sort!(collect(keys(periodic))) == [-1, 1]
  @test maps_agree(periodic[-1], expected)
  @test maps_agree(periodic[1], expected)
  @test harmonics(periodic(t), ω, t) == periodic
end

@testset "non-periodic superoperator phases are rejected" begin
  coefficient = expim(ω * t^2)
  @test_throws ArgumentError harmonics(coefficient * SQA.SuperOperator(dissipator(a)), ω, t)
end

@testset "a zero periodic superoperator keeps its public zero prototype" begin
  periodic = harmonics(zero(SQA.SuperOperator), ω, t)

  @test periodic isa PeriodicGenerator{SQA.SuperOperator}
  @test iszero(periodic)
  @test time_average(periodic) == zero(SQA.SuperOperator)
  @test iszero(
    floquet_expansion(zero(SQA.SuperOperator), ω, t, VanVleck(), 1) |> effective_generator
  )
end

@testset "periodic Lindbladian action agrees with an independent matrix model" begin
  finite_space = NLevelSpace(:finite_validation, 2)
  sigma = Transition(finite_space, :sigma, 1, 2)
  H = sigma + sigma'
  collapse_operator = expim(-ω * t) * sigma
  jump_operator = sigma + expim(ω * t) * sigma'
  rate = 2
  native =
    hamiltonian_action(H) + dissipator(collapse_operator) + rate * dissipator(jump_operator)
  periodic = harmonics(native, ω, t)
  substitutions = Dict(ω => 3.0, t => 0.37)

  H_matrix = tomatrix(H, 2, substitutions)
  C_matrix = tomatrix(collapse_operator, 2, substitutions)
  J_matrix = tomatrix(jump_operator, 2, substitutions)
  rate_value = 2.0
  expected = independent_liouvillian_matrix(H_matrix, C_matrix, J_matrix, rate_value)
  native_matrix = superoperator_matrix(SQA.SuperOperator(native), 2, substitutions)
  periodic_matrix = superoperator_matrix(periodic(t), 2, substitutions)

  @test native_matrix ≈ expected
  @test periodic_matrix ≈ expected
end

@testset "exact drive amplitudes stay exact through harmonics" begin
  h = FockSpace(:c)
  a = Destroy(h, :a)
  drive = (3 // 5) * cos(ω * t) - (4 // 5) * sin(ω * t)
  G = harmonics(hamiltonian_action((2 // 5) * (a + a') * drive), ω, t)
  for (term, c) in G[1]
    left = SQA.left_operator(term)
    expected = isone(left) ? (4 // 25) + (3 // 25) * im : -(4 // 25) - (3 // 25) * im
    @test iszero(SQA.simplify(c - expected))
  end
end

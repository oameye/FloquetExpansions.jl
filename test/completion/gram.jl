using Test
using FloquetExpansions
using LinearAlgebra: eigvals
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables

const SQA = SecondQuantizedAlgebra

function coefficient_matrix_zeros(rows::Int, columns::Int)
  return fill(convert(SQA.CNum, 0), rows, columns)
end

function coefficient_matrix_mul(left, right)
  size(left, 2) == size(right, 1) || throw(DimensionMismatch())
  result = coefficient_matrix_zeros(size(left, 1), size(right, 2))
  for i in axes(result, 1), j in axes(result, 2), k in axes(left, 2)
    result[i, j] = SQA.simplify(result[i, j] + left[i, k] * right[k, j])
  end
  return result
end

function coefficient_matrix_adjoint(matrix)
  result = coefficient_matrix_zeros(size(matrix, 2), size(matrix, 1))
  for i in axes(matrix, 1), j in axes(matrix, 2)
    result[j, i] = conj(matrix[i, j])
  end
  return result
end

function coefficient_matrix_add!(left, right)
  size(left) == size(right) || throw(DimensionMismatch())
  for index in eachindex(left)
    left[index] = SQA.simplify(left[index] + right[index])
  end
  return left
end

function coefficient_matrix_equal(left, right)
  size(left) == size(right) || return false
  return all(iszero(SQA.simplify(left[index] - right[index])) for index in eachindex(left))
end

function gram_coefficient(factorization::GramFactorization, n::Int)
  amplitudes = factorization.amplitudes
  rows = size(first(amplitudes), 1)
  result = coefficient_matrix_zeros(rows, rows)
  for p in 0:n
    p + 1 <= length(amplitudes) || continue
    q = n - p
    q + 1 <= length(amplitudes) || continue
    term = coefficient_matrix_mul(
      amplitudes[p + 1], coefficient_matrix_adjoint(amplitudes[q + 1])
    )
    coefficient_matrix_add!(result, term)
  end
  return result
end

function finite_gram(factorization::GramFactorization)
  amplitudes = factorization.amplitudes
  rows, columns = size(first(amplitudes))
  factor = coefficient_matrix_zeros(rows, columns)
  for amplitude in amplitudes
    coefficient_matrix_add!(factor, amplitude)
  end
  return coefficient_matrix_mul(factor, coefficient_matrix_adjoint(factor))
end

function numeric_matrix(matrix, substitutions)
  result = zeros(ComplexF64, size(matrix))
  for index in eachindex(matrix)
    value = SQA.to_num(matrix[index])
    re = Symbolics.value(Symbolics.substitute(real(value), substitutions))
    im = Symbolics.value(Symbolics.substitute(imag(value), substitutions))
    result[index] = complex(Float64(re), Float64(im))
  end
  return result
end

fock = FockSpace(:gram_completion)
a = Destroy(fock, :a)
@variables ω::Real t::Real γ::Real Δ::Real Ω::Real

@testset "full-rank non-diagonal Gram completion" begin
  frame = DissipativeFrame(a, a^2)
  L = dissipator(a + a^2) + dissipator(a + im * a^2)
  expansion = floquet_expansion(SQA.SuperOperator(L), ω, t, VanVleck(), 2)
  completion = @inferred positive_completion(expansion, Gram(), frame)

  @test dissipative_frame(completion) == frame
  @test factorization(completion) isa GramFactorization
  @test factorization(completion).onsets == [0, 0]
  @test length(factorization(completion).stages) == 1
  @test coefficient_matrix_equal(kossakowski(completion), kossakowski(expansion, frame))
  @test coefficient_matrix_equal(
    finite_gram(factorization(completion)), kossakowski(completion)
  )

  for n in 0:1
    @test coefficient_matrix_equal(
      gram_coefficient(factorization(completion), n), kossakowski_component(completion, n)
    )
    @test effective_component(completion, n) == effective_component(expansion, n)
  end
  @test micromotion(completion) == micromotion(expansion)
  @test SQA.expand_completeness(SQA.SuperOperator(Lindbladian(completion))) ==
    effective_generator(completion)
  @test_throws ArgumentError positive_completion(completion, Gram())
end

@testset "completed Gram accessors return independent containers" begin
  frame = DissipativeFrame(a)
  expansion = floquet_expansion(dissipator(a), ω, t, VanVleck(), 1)
  completion = positive_completion(expansion, Gram(), frame)

  matrix = kossakowski(completion)
  matrix[1, 1] = convert(SQA.CNum, 0)
  @test !iszero(kossakowski(completion)[1, 1])

  completed_channels = channels(completion)
  empty!(completed_channels)
  @test length(channels(completion)) == 1

  gram = factorization(completion)
  gram.amplitudes[1][1, 1] = convert(SQA.CNum, 0)
  @test !iszero(factorization(completion).amplitudes[1][1, 1])
end

canonical_jumps(L) = map(SQA.jump ∘ first, SQA.channels(L))

@testset "automatic frame follows the canonical channel order" begin
  L = dissipator(a + a^2) + dissipator(a + im * a^2)
  expansion = floquet_expansion(L, ω, t, VanVleck(), 1)
  completion = positive_completion(expansion, Gram())

  @test dissipative_frame(completion) == DissipativeFrame(canonical_jumps(L)...)
  @test SQA.expand_completeness(SQA.SuperOperator(Lindbladian(completion))) ==
    effective_generator(completion)

  dependent_L = dissipator(a) + dissipator(a^2) + dissipator(a + a^2)
  dependent = floquet_expansion(dependent_L, ω, t, VanVleck(), 1)
  first_jump, second_jump = canonical_jumps(dependent_L)
  @test dissipative_frame(positive_completion(dependent, Gram())) ==
    DissipativeFrame(first_jump, second_jump)
end

@testset "automatic frame seeds a time-dependent channel with its static harmonics" begin
  n = a' * a
  H = (1//2) * n + (3//10) * n * n + (2//5) * (a + a')
  modulation = 1 + (1//4) * (expim(ω * t) + expim(-ω * t))
  expansion = floquet_expansion(
    hamiltonian_action(H) + (4//5) * dissipator(modulation * a), ω, t, VanVleck(), 3
  )

  @test dissipative_frame(positive_completion(expansion, Gram())) ==
    DissipativeFrame(a, a' * a * a)
end

@testset "symbolic positivity and regularity conditions remain distinct" begin
  expansion = floquet_expansion(γ * dissipator(a), ω, t, VanVleck(), 1)
  completion = positive_completion(expansion, Gram())
  γc = convert(SQA.CNum, γ)

  @test any(iszero(SQA.simplify(value - γc)) for value in positivity_conditions(completion))
  @test any(iszero(SQA.simplify(value - γc)) for value in regularity_conditions(completion))
end

@testset "Cartesian-Pauli driven-qubit completion" begin
  pauli = PauliSpace(:gram_pauli)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σx, σy, σz)
  H = Δ * σz + Ω * cos(ω * t) * σx
  expansion = floquet_expansion(
    hamiltonian_action(H) + dissipator(σx + σy) + dissipator(σy + σz) + dissipator(σz + σx),
    ω,
    t,
    VanVleck(),
    2,
  )
  completion = positive_completion(expansion, Gram(), frame)
  gram = factorization(completion)

  @test size(kossakowski(completion)) == (3, 3)
  @test all(
    coefficient_matrix_equal(
      gram_coefficient(gram, n), kossakowski_component(completion, n)
    ) for n in 0:1
  )
  @test coefficient_matrix_equal(finite_gram(gram), kossakowski(completion))
  @test SQA.expand_completeness(SQA.SuperOperator(Lindbladian(completion))) ==
    effective_generator(completion)

  substitutions = Dict(Δ => 0.7, Ω => 0.4, ω => 10.0)
  eigenvalues = eigvals(numeric_matrix(kossakowski(completion), substitutions))
  @test minimum(real.(eigenvalues)) >= -1.0e-10
end

@testset "bosonic periodically modulated loss" begin
  expansion = floquet_expansion(
    hamiltonian_action(Δ * a' * a) + (2 + cos(ω * t)) * dissipator(a), ω, t, VanVleck(), 2
  )
  completion = positive_completion(expansion, Gram())

  @test dissipative_frame(completion) == DissipativeFrame(a)
  @test kossakowski_component(completion, 0)[1, 1] == 2
  @test iszero(SQA.simplify(kossakowski_component(completion, 1)[1, 1]))
  @test coefficient_matrix_equal(
    finite_gram(factorization(completion)), kossakowski(completion)
  )
end

@testset "retained negative direction is a completion obstruction" begin
  expansion = floquet_expansion(SQA.SuperOperator(-dissipator(a)), ω, t, VanVleck(), 1)
  @test_throws FloquetExpansions.CompletionObstruction positive_completion(
    expansion, Gram(), DissipativeFrame(a)
  )
end

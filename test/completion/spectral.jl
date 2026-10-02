using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using LinearAlgebra: LinearAlgebra

const SQA = SecondQuantizedAlgebra

function spectral_matrix_equal(left, right)
  size(left) == size(right) || return false
  return all(iszero(SQA.simplify(left[index] - right[index])) for index in eachindex(left))
end

function spectral_has_condition(conditions, target)
  coefficient = convert(SQA.CNum, target)
  return any(iszero(SQA.simplify(value - coefficient)) for value in conditions)
end

fock = FockSpace(:spectral_completion)
a = Destroy(fock, :a)
@variables ω::Real t::Real γ::Real Ω::Real

@testset "spectral completion preserves a diagonal rate branch" begin
  frame = DissipativeFrame(a)
  expansion = floquet_expansion(0 * a, ω, t, VanVleck(), 1; channels=(jump(a, γ),))
  completion = @inferred positive_completion(expansion, Spectral(), frame)
  spectral = factorization(completion)

  @test spectral isa SpectralFactorization
  @test spectral.onsets == [0]
  @test spectral.puiseux == [false]
  @test dissipative_frame(completion) == frame
  @test spectral_matrix_equal(kossakowski(completion), kossakowski(expansion, frame))
  @test effective_component(completion, 0) == effective_component(expansion, 0)
  @test micromotion(completion) == micromotion(expansion)
  @test liouvillian(hamiltonian(completion); channels=channels(completion)) ==
    effective_generator(completion)
end

@testset "spectral factorization accessor returns independent containers" begin
  frame = DissipativeFrame(a)
  expansion = floquet_expansion(0 * a, ω, t, VanVleck(), 1; channels=(collapse(a),))
  completion = positive_completion(expansion, Spectral(), frame)

  spectral = factorization(completion)
  spectral.rates[1] = convert(SQA.CNum, 0)
  spectral.vectors[1][1] = convert(SQA.CNum, 0)

  stored = factorization(completion)
  @test !iszero(stored.rates[1])
  @test !iszero(stored.vectors[1][1])
end

@testset "adapted-frame driven-qubit spectral completion" begin
  pauli = PauliSpace(:spectral_pauli)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  expansion = floquet_expansion(
    Ω * cos(ω * t) * σx, ω, t, VanVleck(), 3; channels=(collapse(σz),)
  )
  completion = @inferred positive_completion(expansion, Spectral(), frame)
  spectral = factorization(completion)

  @test spectral isa SpectralFactorization
  @test length(spectral.rates) == 2
  @test length(spectral.vectors) == 2
  @test all(length(vector) == 2 for vector in spectral.vectors)
  for n in 0:2
    @test effective_component(completion, n) == effective_component(expansion, n)
    @test spectral_matrix_equal(
      kossakowski_component(completion, n), kossakowski_component(expansion, frame, n)
    )
  end
  @test micromotion(completion) == micromotion(expansion)
  @test liouvillian(hamiltonian(completion); channels=channels(completion)) ==
    effective_generator(completion)
end

@testset "distinct leading rates generate perturbative branch mixing" begin
  pauli = PauliSpace(:spectral_mixing)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  coherent_rotation = hamiltonian_action(σx)
  diagonal_difference = dissipator(σy) - dissipator(σz)
  leading = dissipator(σz) + 4 * dissipator(σy)
  quadrature = im * (1 // 2) * diagonal_difference
  generator = PeriodicGenerator(
    Dict(
      0 => leading,
      1 => coherent_rotation + quadrature,
      -1 => coherent_rotation - quadrature,
    ),
    ω,
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)
  mixing = kossakowski_component(expansion, frame, 1)

  @test !iszero(SQA.simplify(mixing[1, 2]))
  @test !iszero(SQA.simplify(mixing[2, 1]))

  completion = positive_completion(expansion, Spectral(), frame)
  spectral = factorization(completion)
  @test spectral.onsets == [0, 0]
  @test any(
    !iszero(SQA.simplify(spectral.vectors[column][row])) for
    (column, row) in ((1, 2), (2, 1))
  )
  for n in 0:1
    @test spectral_matrix_equal(
      kossakowski_component(completion, n), kossakowski_component(expansion, frame, n)
    )
  end
end

@testset "degenerate leading sector with retained mixing is rejected" begin
  pauli = PauliSpace(:spectral_degenerate)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  coherent_rotation = hamiltonian_action(σx)
  diagonal_difference = dissipator(σy) - dissipator(σz)
  leading = dissipator(σz) + dissipator(σy)
  quadrature = im * (1 // 2) * diagonal_difference
  generator = PeriodicGenerator(
    Dict(
      0 => leading,
      1 => coherent_rotation + quadrature,
      -1 => coherent_rotation - quadrature,
    ),
    ω,
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)
  mixing = kossakowski_component(expansion, frame, 1)

  @test !iszero(SQA.simplify(mixing[1, 2]))
  @test !iszero(SQA.simplify(mixing[2, 1]))
  @test_throws ArgumentError positive_completion(expansion, Spectral(), frame)
end

@testset "odd positive rate onset remains a valid spectral channel" begin
  pauli = PauliSpace(:spectral_odd_onset)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)

  coherent_rotation = hamiltonian_action(σx)
  cross_dissipator = dissipator(σy + σz) - dissipator(σy) - dissipator(σz)
  quadrature = im * (1 // 2) * cross_dissipator
  generator = PeriodicGenerator(
    Dict(
      0 => dissipator(σz),
      1 => coherent_rotation + quadrature,
      -1 => coherent_rotation - quadrature,
    ),
    ω,
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)
  completion = positive_completion(expansion, Spectral(), frame)
  spectral = factorization(completion)

  @test 1 in spectral.onsets
  odd_branch = findfirst(==(1), spectral.onsets)
  @test odd_branch !== nothing
  @test spectral.puiseux[odd_branch]
  @test spectral_has_condition(positivity_conditions(completion), ω)
  @test spectral_has_condition(regularity_conditions(completion), ω)
  @test liouvillian(hamiltonian(completion); channels=channels(completion)) ==
    effective_generator(completion)
  for n in 0:1
    @test spectral_matrix_equal(
      kossakowski_component(completion, n), kossakowski_component(expansion, frame, n)
    )
  end
end

@testset "Gram and spectral continuations agree through retained qubit order" begin
  pauli = PauliSpace(:spectral_gram_comparison)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  expansion = floquet_expansion(
    Ω * cos(ω * t) * σx, ω, t, VanVleck(), 3; channels=(collapse(σz),)
  )
  gram = positive_completion(expansion, Gram(), frame)
  spectral = positive_completion(expansion, Spectral(), frame)

  for n in 0:2
    retained = kossakowski_component(expansion, frame, n)
    @test spectral_matrix_equal(kossakowski_component(gram, n), retained)
    @test spectral_matrix_equal(kossakowski_component(spectral, n), retained)
  end
  @test effective_component(gram, 2) == effective_component(spectral, 2)
  @test micromotion(gram) == micromotion(spectral)
  @test liouvillian(hamiltonian(gram); channels=channels(gram)) == effective_generator(gram)
  @test liouvillian(hamiltonian(spectral); channels=channels(spectral)) ==
    effective_generator(spectral)
end

@testset "non-diagonal leading spectral frame is rejected" begin
  frame = DissipativeFrame(a, a^2)
  expansion = floquet_expansion(
    0 * a, ω, t, VanVleck(), 1; channels=(collapse(a + a^2), collapse(a + im * a^2))
  )
  @test_throws ArgumentError positive_completion(expansion, Spectral(), frame)
end

function spectral_numeric_matrix(matrix, substitutions)
  return [
    let z = SQA.to_num(entry)
      value(part) = ComplexF64(Symbolics.value(Symbolics.substitute(part, substitutions)))
      value(real(z)) + im * value(imag(z))
    end for entry in matrix
  ]
end

function spectral_numeric_checks(completion, expansion, frame, N, substitutions)
  errors = Float64[]
  for frequency in (8.0, 16.0)
    values = merge(substitutions, Dict(ω => frequency))
    completed = spectral_numeric_matrix(kossakowski(completion), values)
    retained = spectral_numeric_matrix(kossakowski(expansion, frame), values)
    @test minimum(LinearAlgebra.eigvals(LinearAlgebra.Hermitian(completed))) >= -1e-12
    push!(errors, LinearAlgebra.opnorm(completed - retained))
  end
  # The completion changes only terms beyond the retained order N.
  @test log2(errors[1] / errors[2]) >= N + 1 - 0.1
end

@testset "spectral completion of driven qubits beyond the analytic order" begin
  pauli = PauliSpace(:spectral_driven_qubits)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  σm = (1 // 2) * (σx - im * σy)
  σp = (1 // 2) * (σx + im * σy)
  frame = DissipativeFrame(σm, σp, σz)
  @variables E::Real

  # Symbolic decay at order 5: the dark σ₊ branch opens at rate order 4.
  symbolic = floquet_expansion(
    (1 // 2) * σz + E * cos(ω * t) * σx, ω, t, VanVleck(), 5; channels=(jump(σm, γ),)
  )
  completion = positive_completion(symbolic, Spectral(), frame)
  @test factorization(completion).onsets == [0, 4, 2]
  spectral_numeric_checks(completion, symbolic, frame, 4, Dict(E => 0.6, γ => 0.3))

  # Full-rank exact rational rates.
  full_rank = floquet_expansion(
    (1 // 2) * σz + (1 // 2) * cos(ω * t) * σx,
    ω,
    t,
    VanVleck(),
    3;
    channels=(jump(σm, 1 // 5), jump(σp, 1 // 10), jump(σz, 1 // 20)),
  )
  completion = positive_completion(full_rank, Spectral(), frame)
  @test factorization(completion).onsets == [0, 0, 0]
  spectral_numeric_checks(completion, full_rank, frame, 2, Dict{Symbolics.Num,Float64}())

  # Floating-point parameters leave only roundoff outside the frame.
  floating = floquet_expansion(
    (1 // 2) * σz + 0.5 * cos(ω * t) * σx, ω, t, VanVleck(), 5; channels=(jump(σm, 0.2),)
  )
  completion = positive_completion(floating, Spectral(), frame)
  spectral_numeric_checks(completion, floating, frame, 4, Dict{Symbolics.Num,Float64}())
  @test liouvillian(hamiltonian(completion); channels=channels(completion)) ==
    effective_generator(completion)
end

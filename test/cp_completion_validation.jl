using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

function validation_matrix_zeros(rows::Int, columns::Int)
  return fill(convert(SQA.CNum, 0), rows, columns)
end

function validation_matrix_equal(left, right)
  size(left) == size(right) || return false
  return all(iszero(SQA.simplify(left[index] - right[index])) for index in eachindex(left))
end

function validation_matrix_hermitian(matrix)
  size(matrix, 1) == size(matrix, 2) || return false
  return all(
    iszero(SQA.simplify(matrix[row, column] - conj(matrix[column, row]))) for
    row in axes(matrix, 1), column in axes(matrix, 2)
  )
end

function validation_matrix_mul(left, right)
  size(left, 2) == size(right, 1) || throw(DimensionMismatch())
  result = validation_matrix_zeros(size(left, 1), size(right, 2))
  for row in axes(result, 1), column in axes(result, 2), index in axes(left, 2)
    result[row, column] = SQA.simplify(
      result[row, column] + left[row, index] * right[index, column]
    )
  end
  return result
end

function validation_matrix_adjoint(matrix)
  result = validation_matrix_zeros(size(matrix, 2), size(matrix, 1))
  for row in axes(matrix, 1), column in axes(matrix, 2)
    result[column, row] = conj(matrix[row, column])
  end
  return result
end

function validation_matrix_add!(left, right)
  size(left) == size(right) || throw(DimensionMismatch())
  for index in eachindex(left)
    left[index] = SQA.simplify(left[index] + right[index])
  end
  return left
end

function validation_matrix_sub(left, right)
  size(left) == size(right) || throw(DimensionMismatch())
  result = validation_matrix_zeros(size(left)...)
  for index in eachindex(result)
    result[index] = SQA.simplify(left[index] - right[index])
  end
  return result
end

function validation_gram_matrix(factorization::GramFactorization)
  amplitudes = factorization.amplitudes
  rows, columns = size(first(amplitudes))
  factor = validation_matrix_zeros(rows, columns)
  for amplitude in amplitudes
    validation_matrix_add!(factor, amplitude)
  end
  return validation_matrix_mul(factor, validation_matrix_adjoint(factor))
end

function validation_spectral_matrix(factorization::SpectralFactorization)
  q = length(first(factorization.vectors))
  result = validation_matrix_zeros(q, q)
  for branch in eachindex(factorization.rates)
    rate = factorization.rates[branch]
    vector = factorization.vectors[branch]
    for row in 1:q, column in 1:q
      result[row, column] = SQA.simplify(
        result[row, column] + rate * vector[row] * conj(vector[column])
      )
    end
  end
  return result
end

function validation_has_condition(conditions, target; up_to_sign::Bool=false)
  target_coefficient = convert(SQA.CNum, target)
  return any(
    iszero(SQA.simplify(value - target_coefficient)) ||
    (up_to_sign && iszero(SQA.simplify(value + target_coefficient))) for value in conditions
  )
end

function validate_retained_contract(completion, expansion, frame, retained_order)
  @test dissipative_frame(completion) == frame
  for order in 0:retained_order
    @test effective_component(completion, order) == effective_component(expansion, order)
    @test hamiltonian_component(completion, order) == hamiltonian_component(expansion, order)
    @test validation_matrix_equal(
      kossakowski_component(completion, order),
      kossakowski_component(expansion, frame, order),
    )
  end
  @test micromotion(completion) == micromotion(expansion)
  @test validation_matrix_hermitian(kossakowski(completion))
  @test liouvillian(hamiltonian(completion); channels=channels(completion)) ==
    effective_generator(completion)
  return nothing
end

@testset "driven qubit validates Cartesian Gram and adapted spectral frames" begin
  pauli = PauliSpace(:cp_validation_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  bright = σy + σz
  dark = σy - σz
  cartesian = DissipativeFrame(σx, σy, σz)
  adapted = DissipativeFrame(bright, dark)
  @variables ω::Real t::Real Ω::Real

  expansion = floquet_expansion(
    Ω * cos(ω * t) * σx, ω, t, VanVleck(), 3; channels=(collapse(bright),)
  )
  leading_cartesian = kossakowski_component(expansion, cartesian, 0)
  @test !iszero(SQA.simplify(leading_cartesian[2, 3]))
  @test !iszero(SQA.simplify(leading_cartesian[3, 2]))

  gram_cartesian = @inferred positive_completion(expansion, Gram(), cartesian)
  gram_adapted = @inferred positive_completion(expansion, Gram(), adapted)
  spectral_adapted = @inferred positive_completion(expansion, Spectral(), adapted)

  validate_retained_contract(gram_cartesian, expansion, cartesian, 2)
  validate_retained_contract(gram_adapted, expansion, adapted, 2)
  validate_retained_contract(spectral_adapted, expansion, adapted, 2)

  @test validation_matrix_equal(
    validation_gram_matrix(factorization(gram_cartesian)), kossakowski(gram_cartesian)
  )
  @test validation_matrix_equal(
    validation_gram_matrix(factorization(gram_adapted)), kossakowski(gram_adapted)
  )
  @test validation_matrix_equal(
    validation_spectral_matrix(factorization(spectral_adapted)),
    kossakowski(spectral_adapted),
  )

  @test @inferred(dissipative_frame(gram_cartesian)) == cartesian
  gram_channels = @inferred channels(gram_cartesian)
  gram_positivity = @inferred positivity_conditions(gram_cartesian)
  gram_regularity = @inferred regularity_conditions(gram_cartesian)
  gram_factorization = @inferred factorization(gram_cartesian)
  gram_kossakowski = @inferred kossakowski(gram_cartesian)
  @test gram_channels == channels(gram_cartesian)
  @test gram_positivity == positivity_conditions(gram_cartesian)
  @test gram_regularity == regularity_conditions(gram_cartesian)
  @test gram_factorization isa GramFactorization
  @test gram_kossakowski == kossakowski(gram_cartesian)

  @test @inferred(dissipative_frame(spectral_adapted)) == adapted
  spectral_channels = @inferred channels(spectral_adapted)
  spectral_positivity = @inferred positivity_conditions(spectral_adapted)
  spectral_regularity = @inferred regularity_conditions(spectral_adapted)
  spectral_factorization = @inferred factorization(spectral_adapted)
  spectral_kossakowski = @inferred kossakowski(spectral_adapted)
  @test spectral_channels == channels(spectral_adapted)
  @test spectral_positivity == positivity_conditions(spectral_adapted)
  @test spectral_regularity == regularity_conditions(spectral_adapted)
  @test spectral_factorization isa SpectralFactorization
  @test spectral_kossakowski == kossakowski(spectral_adapted)
end

@testset "Gram and spectral completions may differ only beyond retained order" begin
  pauli = PauliSpace(:cp_validation_method_difference)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  @variables ω::Real

  coherent_rotation = hamiltonian_action(σx)
  diagonal_difference = dissipator(σy) - dissipator(σz)
  leading = dissipator(σz) + 4 * dissipator(σy)
  generator = PeriodicGenerator(
    Dict(0 => leading, 1 => coherent_rotation, -1 => diagonal_difference), ω
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)
  gram = positive_completion(expansion, Gram(), frame)
  spectral = positive_completion(expansion, Spectral(), frame)

  validate_retained_contract(gram, expansion, frame, 1)
  validate_retained_contract(spectral, expansion, frame, 1)
  @test !validation_matrix_equal(kossakowski(gram), kossakowski(spectral))
end

@testset "bosonic quadrature-free modulated-loss closure begins beyond retained order" begin
  fock = FockSpace(:cp_validation_modulated_loss)
  a = Destroy(fock, :a)
  frame = DissipativeFrame(a, a^2, a' * a^2)
  @variables ω::Real t::Real K::Real

  H = K * cos(ω * t) * a'^2 * a^2
  expansion = floquet_expansion(
    H,
    ω,
    t,
    VanVleck(),
    2;
    channels=(jump(a, 2 + sin(ω * t)), collapse(a^2)),
  )
  first_order = kossakowski_component(expansion, frame, 1)
  @test any(!iszero(SQA.simplify(first_order[index, 3])) for index in 1:2)
  @test any(!iszero(SQA.simplify(first_order[3, index])) for index in 1:2)

  completion = @inferred positive_completion(expansion, Gram(), frame)
  validate_retained_contract(completion, expansion, frame, 1)
  gram = factorization(completion)
  @test length(gram.amplitudes) == 2
  @test validation_matrix_equal(validation_gram_matrix(gram), kossakowski(completion))

  closure = validation_matrix_sub(kossakowski(completion), kossakowski(expansion, frame))
  second_order_closure = validation_matrix_mul(
    gram.amplitudes[2], validation_matrix_adjoint(gram.amplitudes[2])
  )
  @test validation_matrix_equal(closure, second_order_closure)
  @test !iszero(SQA.simplify(second_order_closure[3, 3]))
end

@testset "Kerr resonator with number-selective loss keeps symbolic positivity" begin
  fock = FockSpace(:cp_validation_kerr)
  a = Destroy(fock, :a)
  number_selective = a' * a^2
  frame = DissipativeFrame(number_selective)
  @variables ω::Real t::Real K::Real γ::Real

  H = K * a'^2 * a^2
  expansion = floquet_expansion(
    H, ω, t, VanVleck(), 1; channels=(jump(number_selective, γ),)
  )

  for method in (Gram(), Spectral())
    completion = @inferred positive_completion(expansion, method, frame)
    validate_retained_contract(completion, expansion, frame, 0)
    @test validation_has_condition(positivity_conditions(completion), γ)
  end
end

@testset "odd positive onset distinguishes Gram amplitudes from spectral rates" begin
  pauli = PauliSpace(:cp_validation_odd_onset)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  @variables ω::Real

  coherent_rotation = hamiltonian_action(σx)
  cross_dissipator = dissipator(σy + σz) - dissipator(σy) - dissipator(σz)
  generator = PeriodicGenerator(
    Dict(0 => dissipator(σz), 1 => coherent_rotation, -1 => cross_dissipator), ω
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)

  @test_throws FractionalJumpOnset positive_completion(expansion, Gram(), frame)
  spectral = positive_completion(expansion, Spectral(), frame)
  spectral_data = factorization(spectral)
  odd_branch = findfirst(==(1), spectral_data.onsets)
  @test odd_branch !== nothing
  @test spectral_data.puiseux[odd_branch]
  validate_retained_contract(spectral, expansion, frame, 1)
end

@testset "retained negative directions use the public obstruction type" begin
  fock = FockSpace(:cp_validation_negative)
  a = Destroy(fock, :a)
  frame = DissipativeFrame(a)
  @variables ω::Real t::Real
  expansion = floquet_expansion(-dissipator(a), ω, t, VanVleck(), 1)

  for method in (Gram(), Spectral())
    error = try
      positive_completion(expansion, method, frame)
      nothing
    catch caught
      caught
    end
    @test error isa CompletionObstruction
    @test error.rate_order == 0
  end
end

@testset "zero diagonal with Hermitian coupling is a PSD obstruction" begin
  fock = FockSpace(:cp_validation_zero_diagonal)
  a = Destroy(fock, :a)
  frame = DissipativeFrame(a, a^2)
  @variables ω::Real t::Real
  cross = dissipator(a + a^2) - dissipator(a) - dissipator(a^2)
  expansion = floquet_expansion(cross, ω, t, VanVleck(), 1)

  error = try
    positive_completion(expansion, Gram(), frame)
    nothing
  catch caught
    caught
  end
  @test error isa CompletionObstruction
  @test error.rate_order == 0
  @test error.reason == :zero_diagonal_coupling
end

@testset "spectral symbolic rank strata separate positivity and regularity" begin
  pauli = PauliSpace(:cp_validation_symbolic_stratum)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  @variables ω::Real γ::Real

  leading = γ * dissipator(σz) + dissipator(σy)
  coherent_rotation = hamiltonian_action(σx)
  diagonal_difference = dissipator(σy) - dissipator(σz)
  generator = PeriodicGenerator(
    Dict(0 => leading, 1 => coherent_rotation, -1 => diagonal_difference), ω
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)
  completion = positive_completion(expansion, Spectral(), frame)

  @test validation_has_condition(positivity_conditions(completion), γ)
  @test validation_has_condition(
    regularity_conditions(completion), γ - 1; up_to_sign=true
  )
  validate_retained_contract(completion, expansion, frame, 1)
end

@testset "Gram completion is covariant under a nonunitary rational frame congruence" begin
  fock = FockSpace(:cp_validation_congruence)
  a = Destroy(fock, :a)
  native = DissipativeFrame(a, a^2)
  g1 = a + (1 // 2) * a^2
  g2 = 2a - a^2
  transformed = DissipativeFrame(g1, g2)
  channel = 2g1 + 3g2
  @variables ω::Real t::Real

  expansion = floquet_expansion(
    0 * a, ω, t, VanVleck(), 1; channels=(collapse(channel),)
  )
  native_completion = positive_completion(expansion, Gram(), native)
  transformed_completion = positive_completion(expansion, Gram(), transformed)

  transform = [
    convert(SQA.CNum, 1) convert(SQA.CNum, 2)
    convert(SQA.CNum, 1 // 2) convert(SQA.CNum, -1)
  ]
  transformed_in_native = validation_matrix_mul(
    transform,
    validation_matrix_mul(
      kossakowski(transformed_completion), validation_matrix_adjoint(transform)
    ),
  )

  @test validation_matrix_equal(kossakowski(native_completion), transformed_in_native)
  @test effective_generator(native_completion) == effective_generator(transformed_completion)
  validate_retained_contract(native_completion, expansion, native, 0)
  validate_retained_contract(transformed_completion, expansion, transformed, 0)
end

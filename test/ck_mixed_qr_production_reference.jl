using Test
using FloquetExpansions
using LinearAlgebra: I

const CKMixedReferenceExact = Complex{Rational{Int}}
const ck_mixed_reference_im = CKMixedReferenceExact(0 // 1, 1 // 1)

function ck_mixed_reference_matrix_unit(dimension::Int, row::Int, column::Int)
  result = zeros(CKMixedReferenceExact, dimension, dimension)
  result[row, column] = one(CKMixedReferenceExact)
  return result
end

function ck_mixed_reference_generator(entries, dimension::Int)
  zero_component = zeros(CKMixedReferenceExact, dimension, dimension)
  terms = Dict{FloquetExpansions.CKBranchVertex{Int},Matrix{CKMixedReferenceExact}}()
  for (vertex, value) in entries
    terms[vertex] = get(terms, vertex, zero_component) + value
  end
  return FloquetExpansions.ck_kernel_generator(terms, zero_component)
end

function ck_mixed_reference_reconstruct(A1, A2, order::Int, dimension::Int)
  identity_component = Matrix{CKMixedReferenceExact}(I, dimension, dimension)
  zero_component = zero(identity_component)
  identity_state = FloquetExpansions.ck_kernel_identity(
    0, identity_component, zero_component
  )
  zero_state = FloquetExpansions.ck_kernel_zero(0, zero_component)
  operations = FloquetExpansions.BlochOrderOperations(
    FloquetExpansions.ck_kernel_product,
    FloquetExpansions.ck_kernel_project_model,
    FloquetExpansions.ck_kernel_solve_complement,
  )
  generators = isnothing(A2) ? [A1] : [A1, A2]
  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    generators, order, identity_state, zero_state, operations
  )
  return FloquetExpansions.evaluate_ck_period_amplitude(
    recurrence.effective, recurrence.wave, order, identity_state, zero_state
  )
end

function ck_mixed_reference_scalarize(polynomial, row::Int, column::Int)
  coefficients = FloquetExpansions.CKOutputKernel{Int,CKMixedReferenceExact}[]
  for kernel in FloquetExpansions.ck_output_kernel(polynomial).coefficients
    terms = Dict{FloquetExpansions.CKOutputKernelKey{Int},CKMixedReferenceExact}()
    for (key, value) in kernel.terms
      coefficient = value[row, column]
      iszero(coefficient) || (terms[key] = coefficient)
    end
    push!(coefficients, FloquetExpansions.CKOutputKernel(terms, zero(CKMixedReferenceExact)))
  end
  return FloquetExpansions.CKOutputPeriodPolynomial(coefficients)
end

function ck_mixed_reference_overlap(left, right)
  return FloquetExpansions.ck_output_channel_pairing(
    left,
    right,
    ck_mixed_reference_im,
    zero(CKMixedReferenceExact),
    (left_value, right_value) -> left_value * conj(right_value),
  )
end

@testset "production endpoint overlap reproduces #197 asymmetric one-output Gram" begin
  dimension = 7
  # Left chronological word: Q_1, R_{a,2}, Q_{-2}, total phase +1.
  # Right chronological word: R_{a,1}, Q_3, total phase +4.
  A1 = ck_mixed_reference_generator(
    [
      FloquetExpansions.ck_jump_vertex(1, 2) =>
        ck_mixed_reference_matrix_unit(dimension, 3, 2),
      FloquetExpansions.ck_jump_vertex(1, 1) =>
        ck_mixed_reference_matrix_unit(dimension, 6, 5),
    ],
    dimension,
  )
  A2 = ck_mixed_reference_generator(
    [
      FloquetExpansions.ck_drift_vertex(1) =>
        ck_mixed_reference_matrix_unit(dimension, 2, 1),
      FloquetExpansions.ck_drift_vertex(-2) =>
        ck_mixed_reference_matrix_unit(dimension, 4, 3),
      FloquetExpansions.ck_drift_vertex(3) =>
        ck_mixed_reference_matrix_unit(dimension, 7, 6),
    ],
    dimension,
  )
  reconstruction = ck_mixed_reference_reconstruct(A1, A2, 5, dimension)
  left = ck_mixed_reference_scalarize(reconstruction.amplitude[6], 4, 1)
  right = ck_mixed_reference_scalarize(reconstruction.amplitude[4], 7, 5)

  forward = ck_mixed_reference_overlap(left, right)
  reverse = ck_mixed_reference_overlap(right, left)
  @test FloquetExpansions.ck_output_pairing_period_terms(forward) ==
    Dict(1 => ck_mixed_reference_im / 6)
  @test FloquetExpansions.ck_output_pairing_period_terms(reverse) ==
    Dict(1 => -ck_mixed_reference_im / 6)
  @test FloquetExpansions.ck_output_pairing_phase_support(forward) == [-3]
  @test FloquetExpansions.ck_output_pairing_phase_support(reverse) == [3]
  @test FloquetExpansions.ck_output_pairing_coefficient(forward, -3, 1) ==
    ck_mixed_reference_im / 6
end

@testset "production endpoint overlap reproduces #197 asymmetric two-output Gram" begin
  dimension = 8
  # Left chronological word: Q_1, R_{a,1}, R_{b,-1}, total phase +1.
  # Right chronological word: R_{a,0}, Q_2, R_{b,-2}, total phase 0.
  A1 = ck_mixed_reference_generator(
    [
      FloquetExpansions.ck_jump_vertex(1, 1) =>
        ck_mixed_reference_matrix_unit(dimension, 3, 2),
      FloquetExpansions.ck_jump_vertex(2, -1) =>
        ck_mixed_reference_matrix_unit(dimension, 4, 3),
      FloquetExpansions.ck_jump_vertex(1, 0) =>
        ck_mixed_reference_matrix_unit(dimension, 6, 5),
      FloquetExpansions.ck_jump_vertex(2, -2) =>
        ck_mixed_reference_matrix_unit(dimension, 8, 7),
    ],
    dimension,
  )
  A2 = ck_mixed_reference_generator(
    [
      FloquetExpansions.ck_drift_vertex(1) =>
        ck_mixed_reference_matrix_unit(dimension, 2, 1),
      FloquetExpansions.ck_drift_vertex(2) =>
        ck_mixed_reference_matrix_unit(dimension, 7, 6),
    ],
    dimension,
  )
  reconstruction = ck_mixed_reference_reconstruct(A1, A2, 4, dimension)
  left = ck_mixed_reference_scalarize(reconstruction.amplitude[5], 4, 1)
  right = ck_mixed_reference_scalarize(reconstruction.amplitude[5], 8, 5)

  forward = ck_mixed_reference_overlap(left, right)
  reverse = ck_mixed_reference_overlap(right, left)
  @test FloquetExpansions.ck_output_pairing_period_terms(forward) ==
    Dict(1 => -(3 // 2) * ck_mixed_reference_im)
  @test FloquetExpansions.ck_output_pairing_period_terms(reverse) ==
    Dict(1 => (3 // 2) * ck_mixed_reference_im)
  @test FloquetExpansions.ck_output_pairing_phase_support(forward) == [1]
  @test FloquetExpansions.ck_output_pairing_phase_support(reverse) == [-1]
end

@testset "production endpoint overlap reproduces #197 pure-jump q=3 limit" begin
  dimension = 8
  A1 = ck_mixed_reference_generator(
    [
      FloquetExpansions.ck_jump_vertex(1, 2) =>
        ck_mixed_reference_matrix_unit(dimension, 2, 1),
      FloquetExpansions.ck_jump_vertex(1, -1) =>
        ck_mixed_reference_matrix_unit(dimension, 3, 2) +
        ck_mixed_reference_matrix_unit(dimension, 4, 3),
      FloquetExpansions.ck_jump_vertex(1, 0) =>
        ck_mixed_reference_matrix_unit(dimension, 6, 5) +
        ck_mixed_reference_matrix_unit(dimension, 7, 6) +
        ck_mixed_reference_matrix_unit(dimension, 8, 7),
    ],
    dimension,
  )
  reconstruction = ck_mixed_reference_reconstruct(A1, nothing, 3, dimension)
  left = ck_mixed_reference_scalarize(reconstruction.amplitude[4], 4, 1)
  right = ck_mixed_reference_scalarize(reconstruction.amplitude[4], 8, 5)

  forward = ck_mixed_reference_overlap(left, right)
  reverse = ck_mixed_reference_overlap(right, left)
  @test FloquetExpansions.ck_output_pairing_period_terms(forward) ==
    Dict(1 => CKMixedReferenceExact(-1 // 2))
  @test FloquetExpansions.ck_output_pairing_period_terms(reverse) ==
    Dict(1 => CKMixedReferenceExact(-1 // 2))
  @test FloquetExpansions.ck_output_pairing_phase_support(forward) == [0]
  @test FloquetExpansions.ck_output_pairing_phase_support(reverse) == [0]
end

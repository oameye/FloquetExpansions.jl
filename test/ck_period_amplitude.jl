using Test
using FloquetExpansions
using LinearAlgebra: I

const CKPeriodExact = Complex{Rational{Int}}
const ck_period_im = CKPeriodExact(0 // 1, 1 // 1)

function ck_period_zero_matrix()
  return zeros(CKPeriodExact, 2, 2)
end

function ck_period_identity_matrix()
  return Matrix{CKPeriodExact}(I, 2, 2)
end

@testset "generic retained-series reconstruction is noncommutative and arbitrary-order" begin
  zero_matrix = ck_period_zero_matrix()
  identity_matrix = ck_period_identity_matrix()
  X1 = CKPeriodExact[0 1; 2 0]
  X2 = CKPeriodExact[1 -1; 0 2]
  B1 = CKPeriodExact[1 2; 0 -1]
  B2 = CKPeriodExact[0 1; -1 1]
  B3 = CKPeriodExact[2 0; 1 -1]

  omega = Matrix{CKPeriodExact}[identity_matrix, X1, X2, zero_matrix]
  omega_inverse = FloquetExpansions.ck_unit_series_inverse(
    omega, 3, identity_matrix, zero_matrix
  )
  identity_series = FloquetExpansions.ck_truncated_series_product(
    omega, omega_inverse, 3, zero_matrix
  )
  @test identity_series == Matrix{CKPeriodExact}[
    identity_matrix, zero_matrix, zero_matrix, zero_matrix
  ]

  generator = Matrix{CKPeriodExact}[zero_matrix, B1, B2, B3]
  slow = FloquetExpansions.ck_zero_constant_series_exponential(
    generator, 3, identity_matrix, zero_matrix
  )
  @test slow[1] == identity_matrix
  @test slow[2] == B1
  @test slow[3] == B2 + (B1 * B1) / 2
  @test slow[4] == B3 + (B1 * B2 + B2 * B1) / 2 + (B1 * B1 * B1) / 6

  dressed = FloquetExpansions.ck_truncated_series_product(
    FloquetExpansions.ck_truncated_series_product(omega, slow, 3, zero_matrix),
    omega_inverse,
    3,
    zero_matrix,
  )
  @test dressed[2] == B1
  @test dressed[3] == B2 + (B1 * B1) / 2 + X1 * B1 - B1 * X1
end

function ck_period_kernel_fixture()
  zero_component = CKPeriodExact(0 // 1, 0 // 1)
  identity_component = CKPeriodExact(1 // 1, 0 // 1)
  components = Dict(
    FloquetExpansions.ck_jump_vertex(1, 0) => identity_component,
    FloquetExpansions.ck_jump_vertex(2, 0) => identity_component,
  )
  A1 = FloquetExpansions.ck_kernel_generator(components, zero_component)
  identity_state = FloquetExpansions.ck_kernel_identity(
    0, identity_component, zero_component
  )
  zero_state = FloquetExpansions.ck_kernel_zero(0, zero_component)
  operations = FloquetExpansions.BlochOrderOperations(
    FloquetExpansions.ck_kernel_product,
    FloquetExpansions.ck_kernel_project_model,
    FloquetExpansions.ck_kernel_solve_complement,
  )
  return (; A1, identity_state, zero_state, operations)
end

function ck_period_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_period_im / mismatch
end

function ck_period_exact_two_output(first_mismatch::Int, second_mismatch::Int)
  result = CKPeriodExact[0, 0, 0]
  if iszero(first_mismatch) && iszero(second_mismatch)
    result[3] = 1 // 2
    return result
  end

  if !iszero(first_mismatch) && iszero(first_mismatch + second_mismatch)
    result[2] += ck_period_inverse_weight(first_mismatch)
  end
  if iszero(first_mismatch) && !iszero(second_mismatch)
    result[2] += ck_period_inverse_weight(second_mismatch)
  elseif !iszero(first_mismatch) && iszero(second_mismatch)
    result[2] -= ck_period_inverse_weight(first_mismatch)
  end
  return result
end

@testset "generic BF endpoint reconstruction equals exact two-output simplex amplitudes" begin
  fixture = ck_period_kernel_fixture()
  recurrence2 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1], 2, fixture.identity_state, fixture.zero_state, fixture.operations
  )
  reconstruction2 = FloquetExpansions.evaluate_ck_period_amplitude(
    recurrence2.effective,
    recurrence2.wave,
    2,
    fixture.identity_state,
    fixture.zero_state,
  )

  @test length(reconstruction2.amplitude) == 3
  two_output = reconstruction2.amplitude[3]
  for first_mismatch in -12:12, second_mismatch in -12:12
    from_bloch = FloquetExpansions.ck_period_ordered_sideband_coefficients(
      two_output,
      [1, 2],
      [-first_mismatch, -second_mismatch];
      inverse_weight=ck_period_inverse_weight,
    )
    @test from_bloch == ck_period_exact_two_output(first_mismatch, second_mismatch)
  end

  @test FloquetExpansions.ck_period_ordered_sideband_coefficients(
    two_output, [1, 2], [0, 0]; inverse_weight=ck_period_inverse_weight
  ) == CKPeriodExact[0, 0, 1 // 2]
  @test FloquetExpansions.ck_period_ordered_sideband_coefficients(
    two_output, [1, 2], [-2, 2]; inverse_weight=ck_period_inverse_weight
  ) == CKPeriodExact[0, ck_period_im / 2, 0]
  @test FloquetExpansions.ck_period_ordered_sideband_coefficients(
    two_output, [1, 2], [0, -3]; inverse_weight=ck_period_inverse_weight
  ) == CKPeriodExact[0, ck_period_im / 3, 0]
  @test FloquetExpansions.ck_period_ordered_sideband_coefficients(
    two_output, [1, 2], [-3, 0]; inverse_weight=ck_period_inverse_weight
  ) == CKPeriodExact[0, -ck_period_im / 3, 0]
  @test all(
    iszero,
    FloquetExpansions.ck_period_ordered_sideband_coefficients(
      two_output, [1, 2], [-2, -3]; inverse_weight=ck_period_inverse_weight
    ),
  )

  recurrence3 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1], 3, fixture.identity_state, fixture.zero_state, fixture.operations
  )
  reconstruction3 = FloquetExpansions.evaluate_ck_period_amplitude(
    recurrence3.effective,
    recurrence3.wave,
    3,
    fixture.identity_state,
    fixture.zero_state,
  )
  @test reconstruction3.amplitude[1:3] == reconstruction2.amplitude
end

using Test
using FloquetExpansions
using LinearAlgebra: I, kron

const CKPeriodRotatingExact = Complex{Rational{Int}}
const ck_period_rotating_im = CKPeriodRotatingExact(0 // 1, 1 // 1)

function ck_period_rotating_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_period_rotating_im / mismatch
end

function ck_period_rotating_exact_overlap(first_mismatch::Int, second_mismatch::Int)
  result = CKPeriodRotatingExact[0, 0, 0]
  if iszero(first_mismatch) && iszero(second_mismatch)
    result[3] = 1 // 2
    return result
  end
  if !iszero(first_mismatch) && iszero(first_mismatch + second_mismatch)
    result[2] += ck_period_rotating_inverse_weight(first_mismatch)
  end
  if iszero(first_mismatch) && !iszero(second_mismatch)
    result[2] += ck_period_rotating_inverse_weight(second_mismatch)
  elseif !iszero(first_mismatch) && iszero(second_mismatch)
    result[2] -= ck_period_rotating_inverse_weight(first_mismatch)
  end
  return result
end

function ck_period_rotating_production_kernel()
  zero_component = CKPeriodRotatingExact(0 // 1, 0 // 1)
  identity_component = CKPeriodRotatingExact(1 // 1, 0 // 1)
  A1 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_jump_vertex(1, 0) => identity_component,
      FloquetExpansions.ck_jump_vertex(2, 0) => identity_component,
    ),
    zero_component,
  )
  identity_state = FloquetExpansions.ck_kernel_identity(
    0, identity_component, zero_component
  )
  zero_state = FloquetExpansions.ck_kernel_zero(0, zero_component)
  operations = FloquetExpansions.BlochOrderOperations(
    FloquetExpansions.ck_kernel_product,
    FloquetExpansions.ck_kernel_project_model,
    FloquetExpansions.ck_kernel_solve_complement,
  )
  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1], 2, identity_state, zero_state, operations
  )
  reconstruction = FloquetExpansions.evaluate_ck_period_amplitude(
    recurrence.effective, recurrence.wave, 2, identity_state, zero_state
  )
  return reconstruction.amplitude[3]
end

function ck_period_rotating_production_overlap(kernel, first_mismatch::Int, second_mismatch::Int)
  return FloquetExpansions.ck_period_ordered_sideband_coefficients(
    kernel,
    [1, 2],
    [-first_mismatch, -second_mismatch];
    inverse_weight=ck_period_rotating_inverse_weight,
  )
end

function ck_period_rotating_fixture()
  sigma_x = CKPeriodRotatingExact[0 1; 1 0]
  sigma_y = CKPeriodRotatingExact[0 -ck_period_rotating_im; ck_period_rotating_im 0]
  sigma_z = CKPeriodRotatingExact[1 0; 0 -1]
  sigma_minus = (1 // 2) * (sigma_x - ck_period_rotating_im * sigma_y)
  sigma_plus = (1 // 2) * (sigma_x + ck_period_rotating_im * sigma_y)
  amplitudes = Dict(
    0 => sigma_z,
    1 => ck_period_rotating_im * sigma_minus,
    -1 => -ck_period_rotating_im * sigma_plus,
  )
  return (; amplitudes, sigma_z)
end

function ck_period_rotating_recycling(
  left::Matrix{CKPeriodRotatingExact}, right::Matrix{CKPeriodRotatingExact}
)
  return kron(conj.(right), left)
end

function ck_period_rotating_hamiltonian_super(H::Matrix{CKPeriodRotatingExact})
  identity_component = Matrix{CKPeriodRotatingExact}(I, 2, 2)
  return -ck_period_rotating_im *
         (kron(identity_component, H) - kron(transpose(H), identity_component))
end

function ck_period_rotating_zero_polynomial(zero_superoperator)
  return Matrix{CKPeriodRotatingExact}[copy(zero_superoperator)]
end

function ck_period_rotating_add_power!(polynomial, power::Int, value, zero_superoperator)
  while length(polynomial) < power + 1
    push!(polynomial, copy(zero_superoperator))
  end
  polynomial[power + 1] += value
  return polynomial
end

function ck_period_rotating_accumulate!(polynomial, coefficients, value, zero_superoperator)
  for power in 0:(length(coefficients) - 1)
    ck_period_rotating_add_power!(
      polynomial, power, coefficients[power + 1] * value, zero_superoperator
    )
  end
  return polynomial
end

function ck_period_rotating_q2_by_phase(amplitudes, zero_superoperator; overlap)
  resolved = Dict{Int,Vector{Matrix{CKPeriodRotatingExact}}}()
  for (first_left, first_left_value) in amplitudes,
    (second_left, second_left_value) in amplitudes,
    (first_right, first_right_value) in amplitudes,
    (second_right, second_right_value) in amplitudes

    coefficients = overlap(first_left - first_right, second_left - second_right)
    all(iszero, coefficients) && continue
    phase = first_left + second_left - first_right - second_right
    polynomial = get!(resolved, phase) do
      return ck_period_rotating_zero_polynomial(zero_superoperator)
    end
    left_value = second_left_value * first_left_value
    right_value = second_right_value * first_right_value
    ck_period_rotating_accumulate!(
      polynomial,
      coefficients,
      ck_period_rotating_recycling(left_value, right_value),
      zero_superoperator,
    )
  end
  return resolved
end

function ck_period_rotating_first_channel(amplitudes, zero_superoperator)
  identity_superoperator = Matrix{CKPeriodRotatingExact}(I, 4, 4)
  recycling_zero = copy(zero_superoperator)
  for value in values(amplitudes)
    recycling_zero += ck_period_rotating_recycling(value, value)
  end
  polynomial = ck_period_rotating_zero_polynomial(zero_superoperator)
  ck_period_rotating_add_power!(
    polynomial, 1, recycling_zero - 2 * identity_superoperator, zero_superoperator
  )
  return polynomial, recycling_zero
end

function ck_period_rotating_full_second_channel(q2_by_phase, recycling_zero, zero_superoperator)
  resolved = Dict(
    phase => [copy(value) for value in polynomial] for (phase, polynomial) in q2_by_phase
  )
  zero_phase = get!(resolved, 0) do
    return ck_period_rotating_zero_polynomial(zero_superoperator)
  end
  identity_superoperator = Matrix{CKPeriodRotatingExact}(I, 4, 4)
  ck_period_rotating_add_power!(
    zero_phase, 2, 2 * identity_superoperator - 2 * recycling_zero, zero_superoperator
  )
  return resolved
end

function ck_period_rotating_polynomial_product(left, right, zero_superoperator)
  result = ck_period_rotating_zero_polynomial(zero_superoperator)
  for left_index in eachindex(left), right_index in eachindex(right)
    ck_period_rotating_add_power!(
      result,
      left_index + right_index - 2,
      left[left_index] * right[right_index],
      zero_superoperator,
    )
  end
  return result
end

function ck_period_rotating_coefficient(polynomial, power::Int, zero_superoperator)
  index = power + 1
  return index <= length(polynomial) ? polynomial[index] : zero_superoperator
end

@testset "production period kernel reproduces the rotating-jump channel/log oracle" begin
  physical_kernel = ck_period_rotating_production_kernel()
  production_overlap = (first, second) ->
    ck_period_rotating_production_overlap(physical_kernel, first, second)
  exact_overlap = ck_period_rotating_exact_overlap

  for first_mismatch in -12:12, second_mismatch in -12:12
    @test production_overlap(first_mismatch, second_mismatch) ==
      exact_overlap(first_mismatch, second_mismatch)
  end

  fixture = ck_period_rotating_fixture()
  zero_superoperator = zeros(CKPeriodRotatingExact, 4, 4)
  q2_production = ck_period_rotating_q2_by_phase(
    fixture.amplitudes, zero_superoperator; overlap=production_overlap
  )
  q2_exact = ck_period_rotating_q2_by_phase(
    fixture.amplitudes, zero_superoperator; overlap=exact_overlap
  )

  @test sort!(collect(keys(q2_production))) == [-2, -1, 0, 1, 2]
  @test keys(q2_production) == keys(q2_exact)
  for phase in keys(q2_exact)
    @test q2_production[phase] == q2_exact[phase]
  end

  first_channel, recycling_zero = ck_period_rotating_first_channel(
    fixture.amplitudes, zero_superoperator
  )
  second_channel = ck_period_rotating_full_second_channel(
    q2_production, recycling_zero, zero_superoperator
  )
  @test iszero(ck_period_rotating_coefficient(second_channel[-2], 1, zero_superoperator))
  @test iszero(ck_period_rotating_coefficient(second_channel[2], 1, zero_superoperator))
  @test !iszero(ck_period_rotating_coefficient(second_channel[-1], 1, zero_superoperator))
  @test !iszero(ck_period_rotating_coefficient(second_channel[1], 1, zero_superoperator))

  first_square = ck_period_rotating_polynomial_product(
    first_channel, first_channel, zero_superoperator
  )
  log_zero_phase = [copy(value) for value in second_channel[0]]
  for power in 0:(length(first_square) - 1)
    ck_period_rotating_add_power!(
      log_zero_phase,
      power,
      -(1 // 2) * ck_period_rotating_coefficient(first_square, power, zero_superoperator),
      zero_superoperator,
    )
  end

  @test iszero(ck_period_rotating_coefficient(log_zero_phase, 2, zero_superoperator))
  expected = ck_period_rotating_hamiltonian_super(-(5 // 4) * fixture.sigma_z)
  @test ck_period_rotating_coefficient(log_zero_phase, 1, zero_superoperator) == expected
end

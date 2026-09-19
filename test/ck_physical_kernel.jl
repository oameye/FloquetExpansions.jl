using Test
using FloquetExpansions
using LinearAlgebra: I

const CKKernelExact = Complex{Rational{Int}}
const ck_kernel_im = CKKernelExact(0 // 1, 1 // 1)

function ck_kernel_fixture()
  zero_component = zeros(CKKernelExact, 2, 2)
  identity_component = Matrix{CKKernelExact}(I, 2, 2)
  jump_components = Dict(
    FloquetExpansions.ck_jump_vertex(1, -1) => CKKernelExact[1 1; 0 -1],
    FloquetExpansions.ck_jump_vertex(1, 0) => CKKernelExact[0 1; 1 1],
    FloquetExpansions.ck_jump_vertex(1, 1) => CKKernelExact[2 -1; 1 0],
  )
  drift_components = Dict(
    FloquetExpansions.ck_drift_vertex(0) => CKKernelExact[1 -1; 2 0],
    FloquetExpansions.ck_drift_vertex(1) => CKKernelExact[0 2; -1 1],
  )
  return (; zero_component, identity_component, jump_components, drift_components)
end

function ck_kernel_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_kernel_im / mismatch
end

function ck_kernel_bosonic_two_output(state, first_sideband::Int, second_sideband::Int)
  result = FloquetExpansions.ck_kernel_ordered_sideband_coefficient(
    state, [first_sideband, second_sideband]; inverse_weight=ck_kernel_inverse_weight
  )
  first_sideband == second_sideband && return result
  return result + FloquetExpansions.ck_kernel_ordered_sideband_coefficient(
    state, [second_sideband, first_sideband]; inverse_weight=ck_kernel_inverse_weight
  )
end

function ck_kernel_commutator_oracle(
  amplitudes, first_sideband::Int, second_sideband::Int, zero_component
)
  result = zero_component
  assignments = first_sideband == second_sideband ?
                ((first_sideband, second_sideband),) :
                ((first_sideband, second_sideband), (second_sideband, first_sideband))

  for (first_vertex, first_value) in amplitudes,
    (second_vertex, second_value) in amplitudes,
    (left_sideband, right_sideband) in assignments

    mismatch = first_vertex.harmonic - left_sideband
    mismatch > 0 || continue
    second_vertex.harmonic + mismatch == right_sideband || continue
    commutator = first_value * second_value - second_value * first_value
    result += (-ck_kernel_im / mismatch) * commutator
  end
  return result
end

@testset "CK physical kernel keeps model/complement and output sidebands distinct" begin
  fixture = ck_kernel_fixture()
  A1 = FloquetExpansions.ck_kernel_generator(
    fixture.jump_components, fixture.zero_component
  )
  model = FloquetExpansions.ck_kernel_project_model(A1)
  solved = FloquetExpansions.ck_kernel_solve_complement(A1)
  solved_twice = FloquetExpansions.ck_kernel_solve_complement(solved)

  @test length(model.terms) == 3
  @test all(key -> key.sector == FloquetExpansions.CKModelSector, keys(model.terms))
  @test length(solved.terms) == 3
  @test all(key -> key.sector == FloquetExpansions.CKComplementSector, keys(solved.terms))
  @test all(key -> key.resolvent_cuts == [1], keys(solved.terms))
  @test all(key -> key.resolvent_cuts == [1, 1], keys(solved_twice.terms))

  for (vertex, value) in fixture.jump_components
    @test FloquetExpansions.ck_kernel_ordered_sideband_coefficient(
      model, [vertex.harmonic]; inverse_weight=ck_kernel_inverse_weight
    ) == value
  end

  for sideband in -8:8
    first_expected = fixture.zero_component
    second_expected = fixture.zero_component
    for (vertex, value) in fixture.jump_components
      mismatch = vertex.harmonic - sideband
      iszero(mismatch) && continue
      first_expected += ck_kernel_inverse_weight(mismatch) * value
      second_expected += ck_kernel_inverse_weight(mismatch)^2 * value
    end
    @test FloquetExpansions.ck_kernel_ordered_sideband_coefficient(
      solved, [sideband]; inverse_weight=ck_kernel_inverse_weight
    ) == first_expected
    @test FloquetExpansions.ck_kernel_ordered_sideband_coefficient(
      solved_twice, [sideband]; inverse_weight=ck_kernel_inverse_weight
    ) == second_expected
  end
end

@testset "generic Bloch recurrence produces the exact finite physical B2 kernel" begin
  fixture = ck_kernel_fixture()
  A1 = FloquetExpansions.ck_kernel_generator(
    fixture.jump_components, fixture.zero_component
  )
  identity_state = FloquetExpansions.ck_kernel_identity(
    0, fixture.identity_component, fixture.zero_component
  )
  zero_state = FloquetExpansions.ck_kernel_zero(0, fixture.zero_component)
  operations = FloquetExpansions.BlochOrderOperations(
    FloquetExpansions.ck_kernel_product,
    FloquetExpansions.ck_kernel_project_model,
    FloquetExpansions.ck_kernel_solve_complement,
  )

  result = @inferred FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1], 2, identity_state, zero_state, operations
  )

  @test result.products == 3
  @test length(result.effective[1].terms) == 3
  @test length(result.wave[1].terms) == 3
  @test length(result.effective[2].terms) == 9
  for key in keys(result.effective[2].terms)
    @test key.sector == FloquetExpansions.CKModelSector
    @test length(key.vertices) == 2
    @test key.resolvent_cuts == [1]
    @test FloquetExpansions.ck_kernel_output_number(key) == 2
  end

  for first_sideband in -8:8, second_sideband in first_sideband:8
    from_kernel = ck_kernel_bosonic_two_output(
      result.effective[2], first_sideband, second_sideband
    )
    from_oracle = ck_kernel_commutator_oracle(
      fixture.jump_components, first_sideband, second_sideband, fixture.zero_component
    )
    @test from_kernel == from_oracle
  end

  far_kernel = ck_kernel_bosonic_two_output(result.effective[2], -32, 32)
  far_oracle = ck_kernel_commutator_oracle(
    fixture.jump_components, -32, 32, fixture.zero_component
  )
  @test far_kernel == far_oracle
  @test !all(iszero, far_kernel)
end

@testset "physical kernel algebra is arbitrary-order and prefix consistent" begin
  fixture = ck_kernel_fixture()
  A1 = FloquetExpansions.ck_kernel_generator(
    fixture.jump_components, fixture.zero_component
  )
  A2 = FloquetExpansions.ck_kernel_generator(
    fixture.drift_components, fixture.zero_component
  )
  identity_state = FloquetExpansions.ck_kernel_identity(
    0, fixture.identity_component, fixture.zero_component
  )
  zero_state = FloquetExpansions.ck_kernel_zero(0, fixture.zero_component)
  operations = FloquetExpansions.BlochOrderOperations(
    FloquetExpansions.ck_kernel_product,
    FloquetExpansions.ck_kernel_project_model,
    FloquetExpansions.ck_kernel_solve_complement,
  )

  order4 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1, A2], 4, identity_state, zero_state, operations
  )
  order5 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1, A2], 5, identity_state, zero_state, operations
  )

  @test order4.effective == order5.effective[1:4]
  @test order4.wave == order5.wave[1:3]
  @test order4.products == 13
  @test order5.products == 19

  for n in eachindex(order5.effective), key in keys(order5.effective[n].terms)
    output_number = FloquetExpansions.ck_kernel_output_number(key)
    drift_number = FloquetExpansions.ck_kernel_drift_number(key)
    @test output_number + 2 * drift_number == n
    @test key.sector == FloquetExpansions.CKModelSector
  end
  for n in eachindex(order5.wave), key in keys(order5.wave[n].terms)
    output_number = FloquetExpansions.ck_kernel_output_number(key)
    drift_number = FloquetExpansions.ck_kernel_drift_number(key)
    @test output_number + 2 * drift_number == n
    @test key.sector == FloquetExpansions.CKComplementSector
  end
end

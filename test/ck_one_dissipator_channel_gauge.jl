using Test
using FloquetExpansions
using LinearAlgebra: I, kron

const CKChannelGaugeExact = Complex{Rational{Int}}
const ck_channel_gauge_im = CKChannelGaugeExact(0 // 1, 1 // 1)

ck_channel_gauge_commutator(left, right) = left * right - right * left

function ck_channel_gauge_fixture()
  sigma_x = CKChannelGaugeExact[0 1; 1 0]
  sigma_z = CKChannelGaugeExact[1 0; 0 -1]

  hamiltonian = Dict(
    0 => sigma_z,
    1 => (1 // 2) * sigma_x,
    -1 => (1 // 2) * sigma_x,
  )
  jumps = Dict(
    0 => sigma_z,
    1 => CKChannelGaugeExact[0 0; 1 0],
  )

  zero_component = zeros(CKChannelGaugeExact, 2, 2)
  identity_component = Matrix{CKChannelGaugeExact}(I, 2, 2)
  A1 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_jump_vertex(1, harmonic) => value for (harmonic, value) in jumps
    ),
    zero_component,
  )
  A2 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_drift_vertex(harmonic) => -ck_channel_gauge_im * value for
      (harmonic, value) in hamiltonian
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
  return (;
    hamiltonian,
    jumps,
    A1,
    A2,
    identity_state,
    zero_state,
    operations,
    zero_component,
  )
end

function ck_channel_gauge_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_channel_gauge_im * (1 // mismatch)
end

function ck_channel_gauge_query(kernel, output_sideband::Int)
  return FloquetExpansions.ck_endpoint_ordered_sideband_coefficient(
    kernel,
    [1],
    [output_sideband];
    inverse_weight=ck_channel_gauge_inverse_weight,
  )
end

function ck_channel_gauge_zero_output_query(kernel)
  return FloquetExpansions.ck_endpoint_ordered_sideband_coefficient(
    kernel, Int[], Int[]; inverse_weight=ck_channel_gauge_inverse_weight
  )
end

function ck_channel_gauge_first_kick(hamiltonian, harmonic, zero_component)
  iszero(harmonic) && return copy(zero_component)
  return ck_channel_gauge_im * (1 // harmonic) *
         get(hamiltonian, harmonic, zero_component)
end

function ck_channel_gauge_second_kick(hamiltonian, harmonic, zero_component)
  iszero(harmonic) && return copy(zero_component)
  Hh = get(hamiltonian, harmonic, zero_component)
  H0 = get(hamiltonian, 0, zero_component)
  result = -ck_channel_gauge_im * (1 // harmonic^2) *
           ck_channel_gauge_commutator(Hh, H0)

  for (inner_harmonic, Hinner) in hamiltonian
    (iszero(inner_harmonic) || inner_harmonic == harmonic) && continue
    Hrest = get(hamiltonian, harmonic - inner_harmonic, zero_component)
    result += -ck_channel_gauge_im * (1 // (2 * harmonic * inner_harmonic)) *
              ck_channel_gauge_commutator(Hinner, Hrest)
  end
  return result
end

function ck_channel_gauge_first_transport(
  hamiltonian, jumps, output_sideband, zero_component
)
  result = copy(zero_component)
  for harmonic in keys(hamiltonian)
    iszero(harmonic) && continue
    jump = get(jumps, output_sideband - harmonic, zero_component)
    kick = ck_channel_gauge_first_kick(hamiltonian, harmonic, zero_component)
    result += ck_channel_gauge_im * ck_channel_gauge_commutator(kick, jump)
  end
  return result
end

function ck_channel_gauge_second_transport(
  hamiltonian, jumps, output_sideband, zero_component
)
  result = copy(zero_component)
  for harmonic in -3:3
    iszero(harmonic) && continue
    jump = get(jumps, output_sideband - harmonic, zero_component)
    kick = ck_channel_gauge_second_kick(hamiltonian, harmonic, zero_component)
    result += ck_channel_gauge_im * ck_channel_gauge_commutator(kick, jump)
  end

  for left_harmonic in keys(hamiltonian), right_harmonic in keys(hamiltonian)
    (iszero(left_harmonic) || iszero(right_harmonic)) && continue
    jump = get(
      jumps, output_sideband - left_harmonic - right_harmonic, zero_component
    )
    left_kick = ck_channel_gauge_first_kick(
      hamiltonian, left_harmonic, zero_component
    )
    right_kick = ck_channel_gauge_first_kick(
      hamiltonian, right_harmonic, zero_component
    )
    nested = ck_channel_gauge_commutator(
      left_kick, ck_channel_gauge_commutator(right_kick, jump)
    )
    result -= (1 // 2) * nested
  end
  return result
end

function ck_channel_gauge_cross_dissipator(left, right)
  identity_component = Matrix{CKChannelGaugeExact}(I, 2, 2)
  norm = adjoint(right) * left
  return kron(conj.(right), left) -
         (1 // 2) *
         (kron(identity_component, norm) + kron(transpose(norm), identity_component))
end

ck_channel_gauge_recycling(left, right) = kron(conj.(right), left)

function ck_channel_gauge_pasted_grade(amplitudes, grade::Int, sidebands; pair)
  zero_superoperator = zeros(CKChannelGaugeExact, 4, 4)
  result = copy(zero_superoperator)
  for left_grade in 0:grade
    right_grade = grade - left_grade
    for sideband in sidebands
      result += pair(
        amplitudes[left_grade + 1][sideband], amplitudes[right_grade + 1][sideband]
      )
    end
  end
  return result
end

function ck_channel_gauge_zero_output_super(amplitude)
  identity_component = Matrix{CKChannelGaugeExact}(I, 2, 2)
  return kron(identity_component, amplitude) + kron(conj.(amplitude), identity_component)
end

function ck_channel_gauge_no_jump_harmonics(jumps, zero_component)
  result = Dict{Int,Matrix{CKChannelGaugeExact}}()
  for (left_harmonic, left) in jumps, (right_harmonic, right) in jumps
    harmonic = left_harmonic - right_harmonic
    value = get(result, harmonic, zero_component) - (1 // 2) * adjoint(right) * left
    if iszero(value)
      haskey(result, harmonic) && delete!(result, harmonic)
    else
      result[harmonic] = value
    end
  end
  return result
end

function ck_channel_gauge_zero_output_B6(fixture, eta::Int)
  no_jump = ck_channel_gauge_no_jump_harmonics(fixture.jumps, fixture.zero_component)
  harmonics = union(keys(fixture.hamiltonian), keys(no_jump))
  components = Dict(
    FloquetExpansions.ck_drift_vertex(harmonic) =>
      -ck_channel_gauge_im * get(fixture.hamiltonian, harmonic, fixture.zero_component) +
      eta * get(no_jump, harmonic, fixture.zero_component) for harmonic in harmonics
  )
  A1_zero = FloquetExpansions.ck_kernel_zero(0, fixture.zero_component)
  A2 = FloquetExpansions.ck_kernel_generator(components, fixture.zero_component)
  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1_zero, A2],
    6,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  canonical = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence.effective,
    recurrence.wave,
    6,
    fixture.identity_state,
    fixture.zero_state,
  )
  return ck_channel_gauge_zero_output_query(canonical.effective[6])
end

function ck_channel_gauge_linear_zero_output_B6(fixture)
  plus_one = ck_channel_gauge_zero_output_B6(fixture, 1)
  minus_one = ck_channel_gauge_zero_output_B6(fixture, -1)
  plus_two = ck_channel_gauge_zero_output_B6(fixture, 2)
  minus_two = ck_channel_gauge_zero_output_B6(fixture, -2)
  return (1 // 12) * (8 * (plus_one - minus_one) - (plus_two - minus_two))
end

function ck_channel_gauge_hamiltonian_super(hamiltonian)
  identity_component = Matrix{CKChannelGaugeExact}(I, 2, 2)
  return -ck_channel_gauge_im *
         (kron(identity_component, hamiltonian) -
          kron(transpose(hamiltonian), identity_component))
end

function ck_channel_gauge_bare_dissipator_harmonic(jumps, harmonic)
  zero_superoperator = zeros(CKChannelGaugeExact, 4, 4)
  result = copy(zero_superoperator)
  for (left_harmonic, left) in jumps, (right_harmonic, right) in jumps
    left_harmonic - right_harmonic == harmonic || continue
    result += ck_channel_gauge_cross_dissipator(left, right)
  end
  return result
end

function ck_channel_gauge_static_similarity(hamiltonian, jumps)
  zero_superoperator = zeros(CKChannelGaugeExact, 4, 4)
  B_R = copy(zero_superoperator)
  for (harmonic, Hh) in hamiltonian
    iszero(harmonic) && continue
    Hmap = ck_channel_gauge_hamiltonian_super(Hh)
    Rminus = ck_channel_gauge_bare_dissipator_harmonic(jumps, -harmonic)
    B_R += (1 // (2 * harmonic^2)) * ck_channel_gauge_commutator(Hmap, Rminus)
  end
  H0 = ck_channel_gauge_hamiltonian_super(hamiltonian[0])
  return B_R, ck_channel_gauge_commutator(B_R, H0)
end

@testset "canonical CK one-dissipator channel matches #104/#122 up to static similarity" begin
  fixture = ck_channel_gauge_fixture()
  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    5,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  canonical_result = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence.effective,
    recurrence.wave,
    5,
    fixture.identity_state,
    fixture.zero_state,
  )
  sidebands = collect(-2:3)

  canonical = [Dict{Int,Matrix{CKChannelGaugeExact}}() for _ in 1:3]
  coherent_cp = [Dict{Int,Matrix{CKChannelGaugeExact}}() for _ in 1:3]
  amplitude_difference = false
  generated_second_order = false

  for sideband in sidebands
    canonical[1][sideband] = ck_channel_gauge_query(
      canonical_result.effective[1], sideband
    )
    canonical[2][sideband] = ck_channel_gauge_query(
      canonical_result.effective[3], sideband
    )
    canonical[3][sideband] = ck_channel_gauge_query(
      canonical_result.effective[5], sideband
    )

    coherent_cp[1][sideband] = get(fixture.jumps, sideband, fixture.zero_component)
    coherent_cp[2][sideband] = ck_channel_gauge_first_transport(
      fixture.hamiltonian, fixture.jumps, sideband, fixture.zero_component
    )
    coherent_cp[3][sideband] = ck_channel_gauge_second_transport(
      fixture.hamiltonian, fixture.jumps, sideband, fixture.zero_component
    )

    @test canonical[1][sideband] == coherent_cp[1][sideband]
    @test canonical[2][sideband] == coherent_cp[2][sideband]
    amplitude_difference |= canonical[3][sideband] != coherent_cp[3][sideband]
    generated_second_order |=
      !haskey(fixture.jumps, sideband) && !iszero(coherent_cp[3][sideband])
  end
  @test amplitude_difference
  @test generated_second_order

  canonical_jump_only_second = ck_channel_gauge_pasted_grade(
    canonical, 2, sidebands; pair=ck_channel_gauge_cross_dissipator
  )
  cp_second = ck_channel_gauge_pasted_grade(
    coherent_cp, 2, sidebands; pair=ck_channel_gauge_cross_dissipator
  )
  canonical_recycling_second = ck_channel_gauge_pasted_grade(
    canonical, 2, sidebands; pair=ck_channel_gauge_recycling
  )
  canonical_zero_output_second = ck_channel_gauge_linear_zero_output_B6(fixture)
  canonical_full_second =
    canonical_recycling_second +
    ck_channel_gauge_zero_output_super(canonical_zero_output_second)

  B_R, expected_similarity = ck_channel_gauge_static_similarity(
    fixture.hamiltonian, fixture.jumps
  )
  @test !iszero(B_R)
  @test !iszero(expected_similarity)

  # The jump-amplitude sector alone supplies only one branch of the map-level gauge.
  @test 2 * (canonical_jump_only_second - cp_second) == expected_similarity
  # Adding the zero-output/no-jump branch gives the complete one-dissipator channel.
  @test canonical_full_second - cp_second == expected_similarity
end

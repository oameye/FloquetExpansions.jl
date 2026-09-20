using Test
using FloquetExpansions
using LinearAlgebra: I, kron

const CKChannelGaugeExact = Complex{Rational{Int}}
const ck_channel_gauge_im = CKChannelGaugeExact(0 // 1, 1 // 1)

ck_channel_gauge_commutator(left, right) = left * right - right * left
ck_channel_gauge_metric_pair(left, right) = adjoint(left) * right
ck_channel_gauge_channel_pair(left, right) = kron(conj.(right), left)

function ck_channel_gauge_fixture()
  sigma_x = CKChannelGaugeExact[0 1; 1 0]
  sigma_z = CKChannelGaugeExact[1 0; 0 -1]
  identity_component = Matrix{CKChannelGaugeExact}(I, 2, 2)
  zero_component = zero(identity_component)

  hamiltonian = Dict(
    0 => sigma_z,
    1 => (1 // 2) * sigma_x,
    -1 => (1 // 2) * sigma_x,
  )
  jumps = Dict(
    0 => sigma_z,
    1 => CKChannelGaugeExact[0 0; 1 0],
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
    identity_component,
    zero_component,
    identity_state,
    zero_state,
    operations,
  )
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

function ck_channel_gauge_amplitudes(fixture, jump_scale::Int, order::Int)
  no_jump = ck_channel_gauge_no_jump_harmonics(fixture.jumps, fixture.zero_component)
  A1 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_jump_vertex(1, harmonic) => jump_scale * value for
      (harmonic, value) in fixture.jumps
    ),
    fixture.zero_component,
  )
  harmonics = union(keys(fixture.hamiltonian), keys(no_jump))
  A2 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_drift_vertex(harmonic) =>
        -ck_channel_gauge_im *
        get(fixture.hamiltonian, harmonic, fixture.zero_component) +
        jump_scale^2 * get(no_jump, harmonic, fixture.zero_component) for
      harmonic in harmonics
    ),
    fixture.zero_component,
  )
  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1, A2],
    order,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  canonical = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence.effective,
    recurrence.wave,
    order,
    fixture.identity_state,
    fixture.zero_state,
  )
  reconstruction = FloquetExpansions.evaluate_ck_period_amplitude(
    canonical.effective,
    canonical.wave,
    order,
    fixture.identity_state,
    fixture.zero_state,
  )
  return [
    FloquetExpansions.ck_output_kernel(amplitude) for amplitude in reconstruction.amplitude
  ]
end

function ck_channel_gauge_connected_order6(fixture, jump_scale::Int)
  order = 6
  amplitudes = ck_channel_gauge_amplitudes(fixture, jump_scale, order)
  metric = FloquetExpansions.ck_output_metric_series(
    amplitudes,
    order,
    ck_channel_gauge_im,
    fixture.zero_component,
    ck_channel_gauge_metric_pair,
  )
  normalization = FloquetExpansions.ck_output_metric_inverse_sqrt_series(
    metric, order, fixture.identity_component
  )
  normalized = FloquetExpansions.ck_output_right_normalize_series(
    amplitudes, normalization, order, fixture.zero_component
  )

  zero_superoperator = zeros(CKChannelGaugeExact, 4, 4)
  channel = FloquetExpansions.ck_output_channel_series(
    normalized,
    order,
    ck_channel_gauge_im,
    zero_superoperator,
    ck_channel_gauge_channel_pair,
  )
  @test all(isempty(channel.coefficients[index].terms) for index in 2:2:6)

  C2 = channel.coefficients[3]
  C4 = channel.coefficients[5]
  C6 = channel.coefficients[7]
  C2C4 = FloquetExpansions.ck_output_pairing_product(C2, C4, zero_superoperator)
  C4C2 = FloquetExpansions.ck_output_pairing_product(C4, C2, zero_superoperator)
  C2C2 = FloquetExpansions.ck_output_pairing_product(C2, C2, zero_superoperator)
  C2C2C2 = FloquetExpansions.ck_output_pairing_product(C2, C2C2, zero_superoperator)

  log6 = FloquetExpansions.CKOutputPairingPolynomial(copy(C6.terms), zero_superoperator)
  FloquetExpansions.ck_output_pairing_add!(
    log6,
    FloquetExpansions.ck_output_pairing_scale(
      -1 // 2,
      FloquetExpansions.ck_output_pairing_add(C2C4, C4C2, zero_superoperator),
      zero_superoperator,
    ),
  )
  FloquetExpansions.ck_output_pairing_add!(
    log6,
    FloquetExpansions.ck_output_pairing_scale(1 // 3, C2C2C2, zero_superoperator),
  )
  return FloquetExpansions.ck_output_pairing_coefficient(log6, 0, 1)
end

function ck_channel_gauge_one_dissipator_ck(fixture)
  value0 = ck_channel_gauge_connected_order6(fixture, 0)
  value1 = ck_channel_gauge_connected_order6(fixture, 1)
  value2 = ck_channel_gauge_connected_order6(fixture, 2)
  value3 = ck_channel_gauge_connected_order6(fixture, 3)
  # At perturbative order six the connected physical generator is an even
  # polynomial in the microscopic jump scale g of degree at most six.  With
  # x=g^2, these are the four interpolation nodes x=0,1,4,9; this combination
  # extracts the coefficient linear in x, i.e. the one-dissipator sector.
  return (-49 // 36) * value0 +
         (3 // 2) * value1 -
         (3 // 20) * value2 +
         (1 // 90) * value3
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

function ck_channel_gauge_cp_second(fixture)
  result = zeros(CKChannelGaugeExact, 4, 4)
  for sideband in -2:3
    leading = get(fixture.jumps, sideband, fixture.zero_component)
    first = ck_channel_gauge_first_transport(
      fixture.hamiltonian, fixture.jumps, sideband, fixture.zero_component
    )
    second = ck_channel_gauge_second_transport(
      fixture.hamiltonian, fixture.jumps, sideband, fixture.zero_component
    )
    result += ck_channel_gauge_cross_dissipator(second, leading)
    result += ck_channel_gauge_cross_dissipator(leading, second)
    result += ck_channel_gauge_cross_dissipator(first, first)
  end
  return result
end

function ck_channel_gauge_hamiltonian_super(hamiltonian)
  identity_component = Matrix{CKChannelGaugeExact}(I, 2, 2)
  return -ck_channel_gauge_im *
         (kron(identity_component, hamiltonian) -
          kron(transpose(hamiltonian), identity_component))
end

function ck_channel_gauge_bare_dissipator_harmonic(jumps, harmonic)
  result = zeros(CKChannelGaugeExact, 4, 4)
  for (left_harmonic, left) in jumps, (right_harmonic, right) in jumps
    left_harmonic - right_harmonic == harmonic || continue
    result += ck_channel_gauge_cross_dissipator(left, right)
  end
  return result
end

function ck_channel_gauge_static_similarity(hamiltonian, jumps)
  B_R = zeros(CKChannelGaugeExact, 4, 4)
  for (harmonic, Hh) in hamiltonian
    iszero(harmonic) && continue
    Hmap = ck_channel_gauge_hamiltonian_super(Hh)
    Rminus = ck_channel_gauge_bare_dissipator_harmonic(jumps, -harmonic)
    B_R += (1 // (2 * harmonic^2)) * ck_channel_gauge_commutator(Hmap, Rminus)
  end
  H0 = ck_channel_gauge_hamiltonian_super(hamiltonian[0])
  return B_R, ck_channel_gauge_commutator(B_R, H0)
end

@testset "fully pasted CK one-dissipator channel matches #104/#122 gauge bridge" begin
  fixture = ck_channel_gauge_fixture()
  canonical_ck_second = ck_channel_gauge_one_dissipator_ck(fixture)
  coherent_cp_second = ck_channel_gauge_cp_second(fixture)
  B_R, expected_similarity = ck_channel_gauge_static_similarity(
    fixture.hamiltonian, fixture.jumps
  )

  @test !iszero(B_R)
  @test !iszero(expected_similarity)
  @test canonical_ck_second != coherent_cp_second
  @test canonical_ck_second - coherent_cp_second == expected_similarity
end

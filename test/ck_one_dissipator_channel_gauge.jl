using Test
using FloquetExpansions
using LinearAlgebra: I, kron

const CKChannelGaugeExact = Complex{Rational{Int}}
const ck_channel_gauge_im = CKChannelGaugeExact(0 // 1, 1 // 1)

function ck_channel_gauge_fixture()
  H0 = CKChannelGaugeExact[2 1 + ck_channel_gauge_im; 1 - ck_channel_gauge_im -1]
  H1 = CKChannelGaugeExact[1 2 - ck_channel_gauge_im; -1 1 + ck_channel_gauge_im]
  H2 = CKChannelGaugeExact[ck_channel_gauge_im 1; 2 -1]
  hamiltonian = Dict(
    0 => H0,
    1 => H1,
    -1 => Matrix(adjoint(H1)),
    2 => H2,
    -2 => Matrix(adjoint(H2)),
  )
  jumps = Dict(
    -1 => CKChannelGaugeExact[1 0; 2 ck_channel_gauge_im],
    0 => CKChannelGaugeExact[0 1; -1 2],
    2 => CKChannelGaugeExact[1 - ck_channel_gauge_im 2; 0 -1],
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
    zero_component,
    identity_component,
    A1,
    A2,
    identity_state,
    zero_state,
    operations,
  )
end

ck_channel_gauge_commutator(left, right) = left * right - right * left

function ck_channel_gauge_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_channel_gauge_im * (1 // mismatch)
end

function ck_channel_gauge_query(kernel, output_sideband::Int)
  return FloquetExpansions.ck_endpoint_ordered_sideband_coefficient(
    kernel, [1], [output_sideband]; inverse_weight=ck_channel_gauge_inverse_weight
  )
end

function ck_channel_gauge_first_kick(hamiltonian, harmonic, zero_component)
  iszero(harmonic) && return copy(zero_component)
  Hh = get(hamiltonian, harmonic, zero_component)
  all(iszero, Hh) && return copy(zero_component)
  return ck_channel_gauge_im * (1 // harmonic) * Hh
end

function ck_channel_gauge_second_kick(hamiltonian, harmonic, zero_component)
  iszero(harmonic) && return copy(zero_component)
  result = copy(zero_component)
  H0 = get(hamiltonian, 0, zero_component)
  Hh = get(hamiltonian, harmonic, zero_component)
  if !all(iszero, Hh)
    result += -ck_channel_gauge_im * (1 // harmonic^2) *
              ck_channel_gauge_commutator(Hh, H0)
  end

  for (inner_harmonic, Hinner) in hamiltonian
    (iszero(inner_harmonic) || inner_harmonic == harmonic) && continue
    Hrest = get(hamiltonian, harmonic - inner_harmonic, zero_component)
    all(iszero, Hrest) && continue
    result += -ck_channel_gauge_im * (1 // (2 * harmonic * inner_harmonic)) *
              ck_channel_gauge_commutator(Hinner, Hrest)
  end
  return result
end

function ck_channel_gauge_second_kick_harmonics(hamiltonian)
  support = Set{Int}()
  harmonics = collect(keys(hamiltonian))
  for harmonic in harmonics
    iszero(harmonic) || push!(support, harmonic)
  end
  for left in harmonics, right in harmonics
    harmonic = left + right
    iszero(harmonic) || push!(support, harmonic)
  end
  return sort!(collect(support))
end

function ck_channel_gauge_first_transport(hamiltonian, jumps, sideband, zero_component)
  result = copy(zero_component)
  for harmonic in keys(hamiltonian)
    iszero(harmonic) && continue
    jump = get(jumps, sideband - harmonic, zero_component)
    all(iszero, jump) && continue
    kick = ck_channel_gauge_first_kick(hamiltonian, harmonic, zero_component)
    result += ck_channel_gauge_im * ck_channel_gauge_commutator(kick, jump)
  end
  return result
end

function ck_channel_gauge_second_transport(hamiltonian, jumps, sideband, zero_component)
  result = copy(zero_component)
  for harmonic in ck_channel_gauge_second_kick_harmonics(hamiltonian)
    jump = get(jumps, sideband - harmonic, zero_component)
    all(iszero, jump) && continue
    kick = ck_channel_gauge_second_kick(hamiltonian, harmonic, zero_component)
    result += ck_channel_gauge_im * ck_channel_gauge_commutator(kick, jump)
  end

  for left_harmonic in keys(hamiltonian), right_harmonic in keys(hamiltonian)
    (iszero(left_harmonic) || iszero(right_harmonic)) && continue
    jump = get(jumps, sideband - left_harmonic - right_harmonic, zero_component)
    all(iszero, jump) && continue
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

function ck_channel_gauge_select_one_output(kernel)
  terms = Dict(
    key => value for (key, value) in kernel.terms if key.output_channels == [1]
  )
  return FloquetExpansions.CKOutputKernel(terms, kernel.zero_component)
end

function ck_channel_gauge_slow_one_output(endpoint_kernel)
  polynomial = FloquetExpansions.ck_period_monomial(endpoint_kernel, 1)
  finite = FloquetExpansions.ck_output_kernel(polynomial)
  coefficients = [ck_channel_gauge_select_one_output(kernel) for kernel in finite.coefficients]
  return FloquetExpansions.CKOutputPeriodPolynomial(coefficients)
end

ck_channel_gauge_channel_pair(left, right) = kron(conj.(right), left)
ck_channel_gauge_metric_pair(left, right) = adjoint(left) * right

function ck_channel_gauge_pasted_pair(left, right, identity_component)
  zero_component = zero(identity_component)
  zero_superoperator = zeros(CKChannelGaugeExact, 4, 4)
  positive = FloquetExpansions.ck_output_channel_pairing(
    left,
    right,
    ck_channel_gauge_im,
    zero_superoperator,
    ck_channel_gauge_channel_pair,
  )
  metric = FloquetExpansions.ck_output_metric_pairing(
    left,
    right,
    ck_channel_gauge_im,
    zero_component,
    ck_channel_gauge_metric_pair,
  )
  positive_matrix = FloquetExpansions.ck_output_pairing_coefficient(positive, 0, 1)
  metric_matrix = FloquetExpansions.ck_output_pairing_coefficient(metric, 0, 1)
  anticommutator =
    (1 // 2) *
    (kron(identity_component, metric_matrix) + kron(transpose(metric_matrix), identity_component))
  return positive_matrix - anticommutator
end

function ck_channel_gauge_pasted_order(amplitudes, pairs, identity_component)
  result = zeros(CKChannelGaugeExact, 4, 4)
  for (left_order, right_order) in pairs
    result += ck_channel_gauge_pasted_pair(
      amplitudes[left_order], amplitudes[right_order], identity_component
    )
  end
  return result
end

function ck_channel_gauge_coherent_sidebands(fixture)
  support = -8:8
  bare = Dict(
    sideband => get(fixture.jumps, sideband, fixture.zero_component) for sideband in support
  )
  first = Dict(
    sideband => ck_channel_gauge_first_transport(
      fixture.hamiltonian, fixture.jumps, sideband, fixture.zero_component
    ) for sideband in support
  )
  second = Dict(
    sideband => ck_channel_gauge_second_transport(
      fixture.hamiltonian, fixture.jumps, sideband, fixture.zero_component
    ) for sideband in support
  )
  return bare, first, second
end

function ck_channel_gauge_coherent_pair(left, right, identity_component)
  positive = zeros(CKChannelGaugeExact, 4, 4)
  metric = zero(identity_component)
  for sideband in keys(left)
    left_value = left[sideband]
    right_value = right[sideband]
    positive += kron(conj.(right_value), left_value)
    metric += adjoint(left_value) * right_value
  end
  anticommutator =
    (1 // 2) *
    (kron(identity_component, metric) + kron(transpose(metric), identity_component))
  return positive - anticommutator
end

function ck_channel_gauge_hamiltonian_super(H, identity_component)
  return -ck_channel_gauge_im *
         (kron(identity_component, H) - kron(transpose(H), identity_component))
end

function ck_channel_gauge_dissipator_harmonic(fixture, harmonic)
  result = zeros(CKChannelGaugeExact, 4, 4)
  for (left_harmonic, left_jump) in fixture.jumps,
    (right_harmonic, right_jump) in fixture.jumps

    left_harmonic - right_harmonic == harmonic || continue
    metric = adjoint(right_jump) * left_jump
    result += kron(conj.(right_jump), left_jump)
    result -= (1 // 2) * kron(fixture.identity_component, metric)
    result -= (1 // 2) * kron(transpose(metric), fixture.identity_component)
  end
  return result
end

function ck_channel_gauge_static_similarity(fixture)
  result = zeros(CKChannelGaugeExact, 4, 4)
  for (harmonic, Hh) in fixture.hamiltonian
    iszero(harmonic) && continue
    Hmap = ck_channel_gauge_hamiltonian_super(Hh, fixture.identity_component)
    Rminus = ck_channel_gauge_dissipator_harmonic(fixture, -harmonic)
    result += (1 // (2 * harmonic^2)) * ck_channel_gauge_commutator(Hmap, Rminus)
  end
  return result
end

@testset "canonical CK one-dissipator gauge bridge holds only after pasting" begin
  fixture = ck_channel_gauge_fixture()
  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    5,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  canonical = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence.effective,
    recurrence.wave,
    5,
    fixture.identity_state,
    fixture.zero_state,
  )

  canonical_amplitudes = Dict(
    0 => ck_channel_gauge_slow_one_output(canonical.effective[1]),
    1 => ck_channel_gauge_slow_one_output(canonical.effective[3]),
    2 => ck_channel_gauge_slow_one_output(canonical.effective[5]),
  )
  coherent0, coherent1, coherent2 = ck_channel_gauge_coherent_sidebands(fixture)

  amplitude_mismatch = false
  for sideband in -7:7
    canonical_second = ck_channel_gauge_query(canonical.effective[5], sideband)
    coherent_second = coherent2[sideband]
    amplitude_mismatch |= canonical_second != coherent_second
  end
  @test amplitude_mismatch

  canonical0 = ck_channel_gauge_pasted_order(
    canonical_amplitudes, [(0, 0)], fixture.identity_component
  )
  canonical1 = ck_channel_gauge_pasted_order(
    canonical_amplitudes, [(0, 1), (1, 0)], fixture.identity_component
  )
  canonical2 = ck_channel_gauge_pasted_order(
    canonical_amplitudes, [(0, 2), (2, 0), (1, 1)], fixture.identity_component
  )

  coherent0_map = ck_channel_gauge_coherent_pair(
    coherent0, coherent0, fixture.identity_component
  )
  coherent1_map =
    ck_channel_gauge_coherent_pair(coherent0, coherent1, fixture.identity_component) +
    ck_channel_gauge_coherent_pair(coherent1, coherent0, fixture.identity_component)
  coherent2_map =
    ck_channel_gauge_coherent_pair(coherent0, coherent2, fixture.identity_component) +
    ck_channel_gauge_coherent_pair(coherent2, coherent0, fixture.identity_component) +
    ck_channel_gauge_coherent_pair(coherent1, coherent1, fixture.identity_component)

  @test canonical0 == coherent0_map
  @test canonical1 == coherent1_map

  static_similarity = ck_channel_gauge_static_similarity(fixture)
  H0map = ck_channel_gauge_hamiltonian_super(
    fixture.hamiltonian[0], fixture.identity_component
  )
  expected_difference = ck_channel_gauge_commutator(static_similarity, H0map)
  @test !iszero(expected_difference)
  @test canonical2 - coherent2_map == expected_difference
end

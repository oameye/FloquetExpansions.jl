using Test
using FloquetExpansions
using LinearAlgebra: I

const CKB5ReductionExact = Complex{Rational{Int}}
const ck_b5_reduction_im = CKB5ReductionExact(0 // 1, 1 // 1)

function ck_b5_reduction_fixture()
  H0 = CKB5ReductionExact[2 1 + ck_b5_reduction_im; 1 - ck_b5_reduction_im -1]
  H1 = CKB5ReductionExact[1 2 - ck_b5_reduction_im; -1 1 + ck_b5_reduction_im]
  H2 = CKB5ReductionExact[ck_b5_reduction_im 1; 2 -1]
  hamiltonian = Dict(
    0 => H0,
    1 => H1,
    -1 => Matrix(adjoint(H1)),
    2 => H2,
    -2 => Matrix(adjoint(H2)),
  )

  jumps = Dict(
    -1 => CKB5ReductionExact[1 0; 2 ck_b5_reduction_im],
    0 => CKB5ReductionExact[0 1; -1 2],
    2 => CKB5ReductionExact[1 - ck_b5_reduction_im 2; 0 -1],
  )

  zero_component = zeros(CKB5ReductionExact, 2, 2)
  identity_component = Matrix{CKB5ReductionExact}(I, 2, 2)
  A1 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_jump_vertex(1, harmonic) => value for (harmonic, value) in jumps
    ),
    zero_component,
  )
  A2 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_drift_vertex(harmonic) => -ck_b5_reduction_im * value for
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
    A1,
    A2,
    identity_state,
    zero_state,
    operations,
  )
end

function ck_b5_reduction_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_b5_reduction_im * (1 // mismatch)
end

function ck_b5_reduction_raw_query(kernel, output_sideband::Int)
  return FloquetExpansions.ck_kernel_ordered_sideband_coefficient(
    kernel, [1], [output_sideband]; inverse_weight=ck_b5_reduction_inverse_weight
  )
end

function ck_b5_reduction_canonical_query(kernel, output_sideband::Int)
  return FloquetExpansions.ck_endpoint_ordered_sideband_coefficient(
    kernel, [1], [output_sideband]; inverse_weight=ck_b5_reduction_inverse_weight
  )
end

ck_b5_reduction_commutator(left, right) = left * right - right * left

function ck_b5_reduction_first_kick(hamiltonian, harmonic, zero_component)
  iszero(harmonic) && return copy(zero_component)
  Hh = get(hamiltonian, harmonic, zero_component)
  all(iszero, Hh) && return copy(zero_component)
  return ck_b5_reduction_im * (1 // harmonic) * Hh
end

function ck_b5_reduction_second_kick(hamiltonian, harmonic, zero_component)
  iszero(harmonic) && return copy(zero_component)
  result = copy(zero_component)
  H0 = get(hamiltonian, 0, zero_component)
  Hh = get(hamiltonian, harmonic, zero_component)
  if !all(iszero, Hh)
    result += -ck_b5_reduction_im * (1 // harmonic^2) *
              ck_b5_reduction_commutator(Hh, H0)
  end

  for (n, Hn) in hamiltonian
    (iszero(n) || n == harmonic) && continue
    Hrest = get(hamiltonian, harmonic - n, zero_component)
    all(iszero, Hrest) && continue
    result += -ck_b5_reduction_im * (1 // (2 * harmonic * n)) *
              ck_b5_reduction_commutator(Hn, Hrest)
  end
  return result
end

function ck_b5_reduction_second_kick_harmonics(hamiltonian)
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

function ck_b5_reduction_first_transport_reference(
  hamiltonian, jumps, output_sideband, zero_component
)
  result = copy(zero_component)
  for harmonic in keys(hamiltonian)
    iszero(harmonic) && continue
    jump = get(jumps, output_sideband - harmonic, zero_component)
    all(iszero, jump) && continue
    kick = ck_b5_reduction_first_kick(hamiltonian, harmonic, zero_component)
    result += ck_b5_reduction_im * ck_b5_reduction_commutator(kick, jump)
  end
  return result
end

function ck_b5_reduction_second_transport_reference(
  hamiltonian, jumps, output_sideband, zero_component
)
  result = copy(zero_component)

  for harmonic in ck_b5_reduction_second_kick_harmonics(hamiltonian)
    jump = get(jumps, output_sideband - harmonic, zero_component)
    all(iszero, jump) && continue
    kick = ck_b5_reduction_second_kick(hamiltonian, harmonic, zero_component)
    result += ck_b5_reduction_im * ck_b5_reduction_commutator(kick, jump)
  end

  for left_harmonic in keys(hamiltonian), right_harmonic in keys(hamiltonian)
    (iszero(left_harmonic) || iszero(right_harmonic)) && continue
    jump = get(
      jumps, output_sideband - left_harmonic - right_harmonic, zero_component
    )
    all(iszero, jump) && continue
    left_kick = ck_b5_reduction_first_kick(
      hamiltonian, left_harmonic, zero_component
    )
    right_kick = ck_b5_reduction_first_kick(
      hamiltonian, right_harmonic, zero_component
    )
    nested = ck_b5_reduction_commutator(
      left_kick, ck_b5_reduction_commutator(right_kick, jump)
    )
    result -= (1 // 2) * nested
  end
  return result
end

@testset "canonical physical B5 one-output sector reduces exactly to #122 second transport" begin
  fixture = ck_b5_reduction_fixture()
  recurrence3 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    3,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  recurrence5 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    5,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  canonical3 = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence3.effective,
    recurrence3.wave,
    3,
    fixture.identity_state,
    fixture.zero_state,
  )
  canonical5 = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence5.effective,
    recurrence5.wave,
    5,
    fixture.identity_state,
    fixture.zero_state,
  )

  sidebands = -7:7

  # B1 and the already-certified first transported amplitude are strict prefixes of the
  # order-five canonical series.
  for sideband in sidebands
    @test ck_b5_reduction_canonical_query(canonical5.effective[1], sideband) ==
      get(fixture.jumps, sideband, fixture.zero_component)
    expected_first = ck_b5_reduction_first_transport_reference(
      fixture.hamiltonian, fixture.jumps, sideband, fixture.zero_component
    )
    @test ck_b5_reduction_canonical_query(canonical5.effective[3], sideband) ==
      expected_first
    @test ck_b5_reduction_canonical_query(canonical3.effective[3], sideband) ==
      expected_first
  end

  # B5^[1]/B1^[1] is two inverse-frequency orders.  After the production
  # P*log(Omega*N)=0 gauge normalization it is exactly #122's second transported jump:
  # i[K_H^(2),L] - 1/2 [K_H^(1),[K_H^(1),L]].
  generated_nonzero = false
  raw_differs = false
  for sideband in sidebands
    expected = ck_b5_reduction_second_transport_reference(
      fixture.hamiltonian, fixture.jumps, sideband, fixture.zero_component
    )
    canonical_value = ck_b5_reduction_canonical_query(
      canonical5.effective[5], sideband
    )
    raw_value = ck_b5_reduction_raw_query(recurrence5.effective[5], sideband)
    @test canonical_value == expected
    raw_differs |= raw_value != canonical_value
    generated_nonzero |= !haskey(fixture.jumps, sideband) && !iszero(expected)
  end
  @test raw_differs
  @test generated_nonzero
end

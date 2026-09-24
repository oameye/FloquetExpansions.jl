using Test
using FloquetExpansions
using LinearAlgebra: I

const CKHDExact = Complex{Rational{Int}}
const ck_hd_im = CKHDExact(0 // 1, 1 // 1)

function ck_hd_fixture()
  H0 = CKHDExact[2 1 + ck_hd_im; 1 - ck_hd_im -1]
  H1 = CKHDExact[1 2 - ck_hd_im; -1 1 + ck_hd_im]
  hamiltonian = Dict(0 => H0, 1 => H1, -1 => Matrix(adjoint(H1)))
  jumps = Dict(
    -1 => CKHDExact[1 0; 2 ck_hd_im],
    0 => CKHDExact[0 1; -1 2],
    1 => CKHDExact[1 - ck_hd_im 2; 0 -1],
  )

  zero_component = zeros(CKHDExact, 2, 2)
  identity_component = Matrix{CKHDExact}(I, 2, 2)
  A1 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_jump_vertex(1, harmonic) => value for (harmonic, value) in jumps
    ),
    zero_component,
  )
  A2 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_drift_vertex(harmonic) => -ck_hd_im * value for
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
  return (; A1, A2, identity_state, zero_state, operations)
end

function ck_hd_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_hd_im * (1 // mismatch)
end

function ck_hd_query(kernel, sidebands::Vector{Int})
  return FloquetExpansions.ck_endpoint_ordered_sideband_coefficient(
    kernel, fill(1, length(sidebands)), sidebands; inverse_weight=ck_hd_inverse_weight
  )
end

function ck_hd_commutator(left, right)
  return FloquetExpansions.ck_endpoint_product(left, right) -
         FloquetExpansions.ck_endpoint_product(right, left)
end

function ck_hd_sideband_samples(outputs::Int)
  outputs == 0 && return [Int[]]
  outputs == 1 && return [[sideband] for sideband in -3:3]
  outputs == 2 && return [[left, right] for left in -2:2 for right in -2:2]
  outputs == 3 &&
    return [[left, middle, right] for left in -1:1 for middle in -1:1 for right in -1:1]
  return [zeros(Int, outputs)]
end

@testset "direct CK Hori-Deprit reproduces the endpoint homological recurrence" begin
  fixture = ck_hd_fixture()
  endpoint_A1 = FloquetExpansions.ck_endpoint_kernel(fixture.A1)
  endpoint_A2 = FloquetExpansions.ck_endpoint_kernel(fixture.A2)
  zero_endpoint = FloquetExpansions.ck_endpoint_kernel(fixture.zero_state)
  hd = FloquetExpansions.evaluate_ck_endpoint_hori_deprit(
    [endpoint_A1, endpoint_A2], 2, zero_endpoint
  )

  expected_G1 = FloquetExpansions.ck_endpoint_solve_homological(endpoint_A1)
  expected_dotG1 = endpoint_A1 - FloquetExpansions.ck_endpoint_project_model(endpoint_A1)
  @test hd.generator[1] == expected_G1
  @test hd.derivative[1] == expected_dotG1

  expected_F2 =
    endpoint_A2 - ck_hd_commutator(expected_G1, endpoint_A1) +
    FloquetExpansions.ck_endpoint_scale(
      ck_hd_commutator(expected_G1, expected_dotG1), 1 // 2
    )
  @test hd.effective[2] == FloquetExpansions.ck_endpoint_project_model(expected_F2)
  @test hd.derivative[2] ==
    expected_F2 - FloquetExpansions.ck_endpoint_project_model(expected_F2)
  @test hd.products > 0
end

@testset "direct CK Hori-Deprit matches canonical Bloch/Feshbach through order five" begin
  fixture = ck_hd_fixture()
  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    5,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  canonical = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence.effective, recurrence.wave, 5, fixture.identity_state, fixture.zero_state
  )
  hd = FloquetExpansions.evaluate_ck_hori_deprit(
    [fixture.A1, fixture.A2], 5, fixture.zero_state
  )

  for order in 1:5
    for key in keys(hd.effective[order].terms)
      outputs = FloquetExpansions.ck_endpoint_output_number(key)
      @test outputs <= order
      @test outputs % 2 == order % 2
    end

    for outputs in 0:order
      outputs % 2 == order % 2 || continue
      for sidebands in ck_hd_sideband_samples(outputs)
        @test ck_hd_query(hd.effective[order], sidebands) ==
          ck_hd_query(canonical.effective[order], sidebands)
      end
    end
  end

  raw_B5 = FloquetExpansions.ck_endpoint_kernel(recurrence.effective[5])
  @test ck_hd_query(raw_B5, [0]) != ck_hd_query(hd.effective[5], [0])
end

@testset "direct CK Hori-Deprit preserves retained prefixes" begin
  fixture = ck_hd_fixture()
  hd4 = FloquetExpansions.evaluate_ck_hori_deprit(
    [fixture.A1, fixture.A2], 4, fixture.zero_state
  )
  hd5 = FloquetExpansions.evaluate_ck_hori_deprit(
    [fixture.A1, fixture.A2], 5, fixture.zero_state
  )

  for order in 1:4
    for outputs in 0:order
      outputs % 2 == order % 2 || continue
      for sidebands in ck_hd_sideband_samples(outputs)
        @test ck_hd_query(hd4.effective[order], sidebands) ==
          ck_hd_query(hd5.effective[order], sidebands)
        @test ck_hd_query(hd4.generator[order], sidebands) ==
          ck_hd_query(hd5.generator[order], sidebands)
      end
    end
  end
  @test hd4.products < hd5.products
end

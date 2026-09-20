using Test
using FloquetExpansions

const CKHDPhysicalExact = Complex{Rational{Int}}
const ck_hd_physical_im = CKHDPhysicalExact(0 // 1, 1 // 1)

function ck_hd_physical_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_hd_physical_im * (1 // mismatch)
end

function ck_hd_physical_query(kernel, sidebands::Vector{Int})
  return FloquetExpansions.ck_endpoint_ordered_sideband_coefficient(
    kernel,
    fill(1, length(sidebands)),
    sidebands;
    inverse_weight=ck_hd_physical_inverse_weight,
  )
end

function ck_hd_physical_b3_fixture()
  H0 = CKHDPhysicalExact[2 1 + ck_hd_physical_im; 1 - ck_hd_physical_im -1]
  H1 = CKHDPhysicalExact[1 2 - ck_hd_physical_im; -1 1 + ck_hd_physical_im]
  H2 = CKHDPhysicalExact[ck_hd_physical_im 1; 2 -1]
  hamiltonian = Dict(
    0 => H0,
    1 => H1,
    -1 => Matrix(adjoint(H1)),
    2 => H2,
    -2 => Matrix(adjoint(H2)),
  )
  jumps = Dict(
    -1 => CKHDPhysicalExact[1 0; 2 ck_hd_physical_im],
    0 => CKHDPhysicalExact[0 1; -1 2],
    2 => CKHDPhysicalExact[1 - ck_hd_physical_im 2; 0 -1],
  )
  zero_component = zeros(CKHDPhysicalExact, 2, 2)
  A1 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_jump_vertex(1, harmonic) => value for (harmonic, value) in jumps
    ),
    zero_component,
  )
  A2 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_drift_vertex(harmonic) => -ck_hd_physical_im * value for
      (harmonic, value) in hamiltonian
    ),
    zero_component,
  )
  zero_state = FloquetExpansions.ck_kernel_zero(0, zero_component)
  return (; hamiltonian, jumps, A1, A2, zero_component, zero_state)
end

function ck_hd_physical_first_transport_reference(
  hamiltonian, jumps, output_sideband, zero_component
)
  result = copy(zero_component)
  for (harmonic, Hh) in hamiltonian
    iszero(harmonic) && continue
    jump = get(jumps, output_sideband - harmonic, zero_component)
    all(iszero, jump) && continue
    kick = ck_hd_physical_im * (1 // harmonic) * Hh
    result += ck_hd_physical_im * (kick * jump - jump * kick)
  end
  return result
end

function ck_hd_physical_rotating_fixture()
  sigma_x = CKHDPhysicalExact[0 1; 1 0]
  sigma_y = CKHDPhysicalExact[0 -ck_hd_physical_im; ck_hd_physical_im 0]
  sigma_z = CKHDPhysicalExact[1 0; 0 -1]
  sigma_minus = (1 // 2) * (sigma_x - ck_hd_physical_im * sigma_y)
  sigma_plus = (1 // 2) * (sigma_x + ck_hd_physical_im * sigma_y)
  amplitudes = Dict(
    0 => sigma_z,
    1 => ck_hd_physical_im * sigma_minus,
    -1 => -ck_hd_physical_im * sigma_plus,
  )
  zero_component = zeros(CKHDPhysicalExact, 2, 2)
  A1 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_jump_vertex(1, harmonic) => value for
      (harmonic, value) in amplitudes
    ),
    zero_component,
  )
  zero_state = FloquetExpansions.ck_kernel_zero(0, zero_component)
  return (; A1, zero_state)
end

@testset "direct CK Hori-Deprit reproduces #263 B3 one-output transport" begin
  fixture = ck_hd_physical_b3_fixture()
  hd = FloquetExpansions.evaluate_ck_hori_deprit(
    [fixture.A1, fixture.A2], 3, fixture.zero_state
  )

  generated_nonzero = false
  for sideband in -4:4
    expected = ck_hd_physical_first_transport_reference(
      fixture.hamiltonian, fixture.jumps, sideband, fixture.zero_component
    )
    actual = ck_hd_physical_query(hd.effective[3], [sideband])
    @test actual == expected
    generated_nonzero |= !haskey(fixture.jumps, sideband) && !iszero(expected)
  end
  @test generated_nonzero
end

@testset "direct CK Hori-Deprit reproduces rotating-jump B2 two-output kernel" begin
  fixture = ck_hd_physical_rotating_fixture()
  endpoint_A1 = FloquetExpansions.ck_endpoint_kernel(fixture.A1)
  zero_endpoint = FloquetExpansions.ck_endpoint_kernel(fixture.zero_state)
  hd = FloquetExpansions.evaluate_ck_endpoint_hori_deprit([endpoint_A1], 2, zero_endpoint)
  expected_B2 = FloquetExpansions.ck_endpoint_project_model(
    FloquetExpansions.ck_endpoint_product(
      endpoint_A1, FloquetExpansions.ck_endpoint_solve_homological(endpoint_A1)
    )
  )

  nonzero = false
  for left_sideband in -2:2, right_sideband in -2:2
    sidebands = [left_sideband, right_sideband]
    expected = ck_hd_physical_query(expected_B2, sidebands)
    actual = ck_hd_physical_query(hd.effective[2], sidebands)
    @test actual == expected
    nonzero |= !iszero(expected)
  end
  @test nonzero
  @test all(
    FloquetExpansions.ck_endpoint_output_number(key) == 2 for key in keys(hd.effective[2].terms)
  )
end

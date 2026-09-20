using Test
using FloquetExpansions
using LinearAlgebra: I

const CKOrder6Exact = Complex{Rational{Int}}
const ck_order6_im = CKOrder6Exact(0 // 1, 1 // 1)

function ck_order6_fixture()
  sigma_x = CKOrder6Exact[0 1; 1 0]
  sigma_z = CKOrder6Exact[1 0; 0 -1]
  identity_component = Matrix{CKOrder6Exact}(I, 2, 2)
  zero_component = zero(identity_component)
  hamiltonian = Dict(0 => sigma_z, 1 => (1 // 2) * sigma_x, -1 => (1 // 2) * sigma_x)
  jumps = Dict(0 => sigma_z, 1 => CKOrder6Exact[0 0; 1 0])

  no_jump = Dict{Int,Matrix{CKOrder6Exact}}()
  for (left_harmonic, left) in jumps, (right_harmonic, right) in jumps
    harmonic = left_harmonic - right_harmonic
    value = get(no_jump, harmonic, zero_component) - (1 // 2) * adjoint(right) * left
    if iszero(value)
      haskey(no_jump, harmonic) && delete!(no_jump, harmonic)
    else
      no_jump[harmonic] = value
    end
  end

  A1 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_jump_vertex(1, harmonic) => value for (harmonic, value) in jumps
    ),
    zero_component,
  )
  harmonics = union(keys(hamiltonian), keys(no_jump))
  A2 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_drift_vertex(harmonic) =>
        -ck_order6_im * get(hamiltonian, harmonic, zero_component) +
        get(no_jump, harmonic, zero_component) for harmonic in harmonics
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

function ck_order6_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_order6_im * (1 // mismatch)
end

function ck_order6_endpoint_query(kernel, sidebands::Vector{Int})
  return FloquetExpansions.ck_endpoint_ordered_sideband_coefficient(
    kernel, fill(1, length(sidebands)), sidebands; inverse_weight=ck_order6_inverse_weight
  )
end

function ck_order6_output_query(polynomial, sidebands::Vector{Int})
  return FloquetExpansions.ck_output_ordered_sideband_coefficients(
    polynomial,
    fill(1, length(sidebands)),
    sidebands;
    inverse_weight=ck_order6_inverse_weight,
  )
end

function ck_order6_sideband_samples(outputs::Int)
  outputs == 0 && return [Int[]]
  outputs == 1 && return [[sideband] for sideband in -3:3]
  outputs == 2 && return [[left, right] for left in -2:2 for right in -2:2]
  outputs == 3 &&
    return [[left, middle, right] for left in -1:1 for middle in -1:1 for right in -1:1]
  return [zeros(Int, outputs)]
end

@testset "one-dissipator CK order-six physical gauge boundary" begin
  fixture = ck_order6_fixture()
  order = 6
  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    order,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  canonical = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence.effective, recurrence.wave, order, fixture.identity_state, fixture.zero_state
  )
  hd = FloquetExpansions.evaluate_ck_hori_deprit(
    [fixture.A1, fixture.A2], order, fixture.zero_state
  )

  # #281 certifies canonical Bloch/Feshbach = direct Hori-Deprit only through
  # order five.  The connected channel theorem consumes the order-six slow
  # coefficient as well as the first five embedding coefficients.  Compare
  # physical ordered-sideband semantics rather than private kernel encoding.
  for perturbative_order in 1:5
    for outputs in 0:perturbative_order
      outputs % 2 == perturbative_order % 2 || continue
      for sidebands in ck_order6_sideband_samples(outputs)
        @test ck_order6_endpoint_query(
          canonical.log_embedding[perturbative_order], sidebands
        ) == ck_order6_endpoint_query(hd.generator[perturbative_order], sidebands)
      end
    end
  end

  for perturbative_order in 1:6
    for outputs in 0:perturbative_order
      outputs % 2 == perturbative_order % 2 || continue
      for sidebands in ck_order6_sideband_samples(outputs)
        @test ck_order6_endpoint_query(
          canonical.effective[perturbative_order], sidebands
        ) == ck_order6_endpoint_query(hd.effective[perturbative_order], sidebands)
      end
    end
  end

  identity_endpoint = FloquetExpansions.ck_endpoint_kernel(fixture.identity_state)
  zero_endpoint = FloquetExpansions.ck_endpoint_kernel(fixture.zero_state)
  raw_period = FloquetExpansions.evaluate_ck_period_amplitude(
    [FloquetExpansions.ck_endpoint_kernel(value) for value in recurrence.effective],
    [FloquetExpansions.ck_endpoint_kernel(value) for value in recurrence.wave],
    order,
    identity_endpoint,
    zero_endpoint,
  )
  canonical_period = FloquetExpansions.evaluate_ck_period_amplitude(
    canonical.effective, canonical.wave, order, identity_endpoint, zero_endpoint
  )
  raw_output = [FloquetExpansions.ck_output_kernel(value) for value in raw_period.amplitude]
  canonical_output = [
    FloquetExpansions.ck_output_kernel(value) for value in canonical_period.amplitude
  ]

  # The static canonical model-space similarity changes the internal endpoint
  # representation but not the physical one-period amplitude when the wave
  # and slow generator are transformed together.
  for perturbative_order in 0:order
    for outputs in 0:perturbative_order
      outputs % 2 == perturbative_order % 2 || continue
      for sidebands in ck_order6_sideband_samples(outputs)
        @test ck_order6_output_query(raw_output[perturbative_order + 1], sidebands) ==
          ck_order6_output_query(canonical_output[perturbative_order + 1], sidebands)
      end
    end
  end
end

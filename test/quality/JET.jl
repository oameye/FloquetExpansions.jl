using Test
using FloquetExpansions
using JET: JET
using Symbolics: @variables

@testset "JET" begin
  JET.test_package(FloquetExpansions; target_modules=(FloquetExpansions,))
end

@testset "expansion algorithm optimizer stability" begin
  pauli = PauliSpace(:jet_expansion_algorithm)
  σx = Pauli(pauli, :sigma, 1)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω_bf::Real
  generator = PeriodicGenerator(Dict(0 => 1 * σz, 1 => 1 * σx, -1 => 1 * σx), ω_bf)

  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    generator, VanVleck(), 3
  )
  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    generator, VanVleck(; algorithm=BlochFeshbach()), 3
  )

  liouvillian_generator = PeriodicGenerator(
    Dict(0 => hamiltonian_action(σz), 1 => dissipator(σx), -1 => hamiltonian_action(σx)),
    ω_bf,
  )
  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    liouvillian_generator, VanVleck(; algorithm=BlochFeshbach()), 2
  )
end

@testset "CP algorithm policy optimizer stability" begin
  JET.@test_opt target_modules=(FloquetExpansions,) HoriDeprit(;
    complete_positive=Val(false)
  )
  JET.@test_opt target_modules=(FloquetExpansions,) BlochFeshbach(;
    complete_positive=Val(false)
  )
  JET.@test_opt target_modules=(FloquetExpansions,) HoriDeprit(;
    complete_positive=Val(true)
  )
  JET.@test_opt target_modules=(FloquetExpansions,) BlochFeshbach(;
    complete_positive=Val(true)
  )

  pauli = PauliSpace(:jet_cp_algorithm_policy)
  σx = Pauli(pauli, :sigma, 1)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω_cp_algorithm_policy::Real
  hamiltonian_generator = PeriodicGenerator(
    Dict(0 => 1 * σz, 1 => 1 * σx, -1 => 1 * σx), ω_cp_algorithm_policy
  )
  liouvillian_generator = PeriodicGenerator(
    Dict(0 => hamiltonian_action(σz) + dissipator(σx)), ω_cp_algorithm_policy
  )

  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    hamiltonian_generator,
    VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(true))),
    2,
  )
  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    hamiltonian_generator,
    VanVleck(; algorithm=BlochFeshbach(; complete_positive=Val(true))),
    2,
  )
  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    liouvillian_generator,
    VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))),
    2,
  )
  JET.@test_opt target_modules=(FloquetExpansions,) floquet_expansion(
    liouvillian_generator,
    VanVleck(; algorithm=BlochFeshbach(; complete_positive=Val(false))),
    2,
  )
end

@testset "Bloch order recurrence optimizer stability" begin
  operations = FloquetExpansions.BlochOrderOperations(*, identity, identity)
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.evaluate_bloch_order_recurrence(
    [1.0, 0.25], 5, 1.0, 0.0, operations
  )
end

@testset "physical CK kernel optimizer stability" begin
  jumps = Dict(
    FloquetExpansions.ck_jump_vertex(1, -1) => 1.0,
    FloquetExpansions.ck_jump_vertex(1, 0) => 2.0,
    FloquetExpansions.ck_jump_vertex(1, 1) => 3.0,
  )
  drifts = Dict(FloquetExpansions.ck_drift_vertex(0) => -0.5)
  A1 = FloquetExpansions.ck_kernel_generator(jumps, 0.0)
  A2 = FloquetExpansions.ck_kernel_generator(drifts, 0.0)
  identity_state = FloquetExpansions.ck_kernel_identity(0, 1.0, 0.0)
  zero_state = FloquetExpansions.ck_kernel_zero(0, 0.0)
  operations = FloquetExpansions.BlochOrderOperations(
    FloquetExpansions.ck_kernel_product,
    FloquetExpansions.ck_kernel_project_model,
    FloquetExpansions.ck_kernel_solve_complement,
  )

  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1, A2], 4, identity_state, zero_state, operations
  )

  multichannel = FloquetExpansions.ck_kernel_project_model(
    FloquetExpansions.ck_kernel_generator(
      Dict(
        FloquetExpansions.ck_jump_vertex(1, 0) => 2.0,
        FloquetExpansions.ck_jump_vertex(2, 0) => 3.0,
      ),
      0.0,
    ),
  )
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.ck_kernel_ordered_sideband_coefficient(
    multichannel, [2], [0]; inverse_weight=inv
  )

  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1, A2], 3, identity_state, zero_state, operations
  )
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.evaluate_ck_period_amplitude(
    recurrence.effective, recurrence.wave, 3, identity_state, zero_state
  )
  reconstruction = FloquetExpansions.evaluate_ck_period_amplitude(
    recurrence.effective, recurrence.wave, 3, identity_state, zero_state
  )
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.ck_period_ordered_sideband_coefficients(
    reconstruction.amplitude[3], [1, 1], [-1, 1]; inverse_weight=inv
  )

  recurrence5 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1, A2], 5, identity_state, zero_state, operations
  )
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence5.effective, recurrence5.wave, 5, identity_state, zero_state
  )

  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.evaluate_ck_hori_deprit(
    [A1, A2], 3, zero_state
  )
end

@testset "physical CK normalization optimizer stability" begin
  metric = FloquetExpansions.CKOutputPairingSeries([
    FloquetExpansions.CKOutputPairingPolynomial(Dict(0 => 1.0), 0.0),
    FloquetExpansions.CKOutputPairingPolynomial(Dict(1 => 0.25), 0.0),
  ])
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.ck_output_metric_inverse_sqrt_series(
    metric, 4, 1.0
  )

  normalization = FloquetExpansions.ck_output_metric_inverse_sqrt_series(metric, 4, 1.0)
  vacuum = FloquetExpansions.CKOutputKernelKey(
    Int[],
    FloquetExpansions.CKOutputBlock[],
    FloquetExpansions.CKOutputConstraint{Int}[],
    FloquetExpansions.CKOutputConstraint{Int}[],
  )
  amplitude = FloquetExpansions.CKOutputPeriodPolynomial([
    FloquetExpansions.CKOutputKernel(Dict(vacuum => 1.0), 0.0)
  ])
  JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.ck_output_right_normalize_series(
    [amplitude], normalization, 4, 0.0
  )
end

@testset "completion optimizer stability" begin
  fock = FockSpace(:jet_completion_fock)
  a = Destroy(fock, :a)
  @variables ω::Real t::Real
  gram_frame = DissipativeFrame(a, a^2)
  gram_generator = liouvillian(0 * a; channels=(collapse(a + a^2), collapse(a + im * a^2)))
  gram_expansion = floquet_expansion(gram_generator, ω, t, VanVleck(), 1)

  JET.@test_opt target_modules=(FloquetExpansions,) positive_completion(
    gram_expansion, Gram(), gram_frame
  )

  pauli = PauliSpace(:jet_completion_pauli)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  recursive_frame = DissipativeFrame(σx, σy, σz)
  recursive_expansion = floquet_expansion(
    cos(ω * t) * σx, ω, t, VanVleck(), 3; channels=(collapse(σz),)
  )

  JET.@test_opt target_modules=(FloquetExpansions,) positive_completion(
    recursive_expansion, Gram(), recursive_frame
  )

  spectral_frame = DissipativeFrame(σz, σy)
  JET.@test_opt target_modules=(FloquetExpansions,) positive_completion(
    recursive_expansion, Spectral(), spectral_frame
  )

  completion = positive_completion(gram_expansion, Gram(), gram_frame)
  JET.@test_opt target_modules=(FloquetExpansions,) hamiltonian(completion)
  JET.@test_opt target_modules=(FloquetExpansions,) kossakowski_component(completion, 0)
end

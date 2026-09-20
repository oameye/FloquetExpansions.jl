using Test
using FloquetExpansions
using LinearAlgebra: I

const CKB3ReductionExact = Complex{Rational{Int}}
const ck_b3_reduction_im = CKB3ReductionExact(0 // 1, 1 // 1)

function ck_b3_reduction_fixture()
  H0 = CKB3ReductionExact[2 1 + ck_b3_reduction_im; 1 - ck_b3_reduction_im -1]
  H1 = CKB3ReductionExact[1 2 - ck_b3_reduction_im; -1 1 + ck_b3_reduction_im]
  H2 = CKB3ReductionExact[ck_b3_reduction_im 1; 2 -1]
  hamiltonian = Dict(
    0 => H0,
    1 => H1,
    -1 => adjoint(H1),
    2 => H2,
    -2 => adjoint(H2),
  )

  jumps = Dict(
    -1 => CKB3ReductionExact[1 0; 2 ck_b3_reduction_im],
    0 => CKB3ReductionExact[0 1; -1 2],
    2 => CKB3ReductionExact[1 - ck_b3_reduction_im 2; 0 -1],
  )

  zero_component = zeros(CKB3ReductionExact, 2, 2)
  identity_component = Matrix{CKB3ReductionExact}(I, 2, 2)
  A1 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_jump_vertex(1, harmonic) => value for
      (harmonic, value) in jumps
    ),
    zero_component,
  )
  A2 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_drift_vertex(harmonic) => -ck_b3_reduction_im * value for
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

function ck_b3_reduction_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_b3_reduction_im * (1 // mismatch)
end

function ck_b3_reduction_query(kernel, output_sideband::Int)
  return FloquetExpansions.ck_kernel_ordered_sideband_coefficient(
    kernel,
    [1],
    [output_sideband];
    inverse_weight=ck_b3_reduction_inverse_weight,
  )
end

function ck_b3_reduction_vacuum_query(kernel)
  return FloquetExpansions.ck_kernel_ordered_sideband_coefficient(
    kernel,
    Int[],
    Int[];
    inverse_weight=ck_b3_reduction_inverse_weight,
  )
end

function ck_b3_reduction_transport_reference(hamiltonian, jumps, output_sideband, zero_component)
  result = copy(zero_component)
  for (harmonic, Hh) in hamiltonian
    iszero(harmonic) && continue
    jump = get(jumps, output_sideband - harmonic, zero_component)
    all(iszero, jump) && continue

    # #122 convention: K_H^(1)[h] = i H_h/h and
    # A_L^(1) = i[K_H^(1), L].
    kick = ck_b3_reduction_im * (1 // harmonic) * Hh
    result += ck_b3_reduction_im * (kick * jump - jump * kick)
  end
  return result
end

@testset "physical B3 one-output sector reduces exactly to #122 amplitude transport" begin
  fixture = ck_b3_reduction_fixture()
  recurrence3 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    3,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  recurrence4 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    4,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )

  B1 = recurrence3.effective[1]
  B2 = recurrence3.effective[2]
  B3 = recurrence3.effective[3]
  X1 = recurrence3.wave[1]
  X2 = recurrence3.wave[2]

  # Leading one-output model amplitudes are exactly the bare physical jump rows.
  for sideband in -4:4
    @test ck_b3_reduction_query(B1, sideband) ==
      get(fixture.jumps, sideband, fixture.zero_component)
  end

  # In the strict one-dissipator/Born slice A2=-iH, the zero-output slow block is -iH0.
  @test ck_b3_reduction_vacuum_query(B2) ==
    -ck_b3_reduction_im * fixture.hamiltonian[0]

  # The two BF folds cannot return the one-output complement branch to the model space.
  fold_x1_b2 = FloquetExpansions.ck_kernel_project_model(
    FloquetExpansions.ck_kernel_product(X1, B2)
  )
  fold_x2_b1 = FloquetExpansions.ck_kernel_project_model(
    FloquetExpansions.ck_kernel_product(X2, B1)
  )
  for sideband in -4:4
    @test iszero(ck_b3_reduction_query(fold_x1_b2, sideband))
    @test iszero(ck_b3_reduction_query(fold_x2_b1, sideband))
  end

  # B3^[1]/B1^[1] is one inverse-frequency order.  Sideband by sideband it is precisely
  # the first #122 transported-jump correction i[K_H^(1),L].  The range includes generated
  # sidebands outside the bare jump support.
  for sideband in -4:4
    expected = ck_b3_reduction_transport_reference(
      fixture.hamiltonian, fixture.jumps, sideband, fixture.zero_component
    )
    @test ck_b3_reduction_query(B3, sideband) == expected
  end
  @test !iszero(ck_b3_reduction_query(B3, 4))
  @test !haskey(fixture.jumps, 4)

  # H0 cannot enter K_H^(1), and the production recurrence realizes the same cancellation.
  hamiltonian_without_zero = copy(fixture.hamiltonian)
  delete!(hamiltonian_without_zero, 0)
  for sideband in -4:4
    @test ck_b3_reduction_transport_reference(
      fixture.hamiltonian, fixture.jumps, sideband, fixture.zero_component
    ) == ck_b3_reduction_transport_reference(
      hamiltonian_without_zero, fixture.jumps, sideband, fixture.zero_component
    )
  end

  # Requesting one higher perturbative order must not change the already-retained B3 data.
  @test recurrence4.effective[3] == B3
  @test recurrence4.wave[1] == recurrence3.wave[1]
  @test recurrence4.wave[2] == recurrence3.wave[2]

  grades = sort!(unique(FloquetExpansions.ck_kernel_output_number(key) for key in keys(B3.terms)))
  @test 1 in grades
  @test all(grade in (1, 3) for grade in grades)
end

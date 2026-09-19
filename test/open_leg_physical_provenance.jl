using Test
using LinearAlgebra: I, kron

include(joinpath(@__DIR__, "helpers", "open_leg_reference.jl"))
include(joinpath(@__DIR__, "helpers", "open_leg_evaluation_plan.jl"))
include(joinpath(@__DIR__, "helpers", "stinespring_normalization_reference.jl"))
include(joinpath(@__DIR__, "helpers", "open_leg_physical_provenance.jl"))

const PhysicalExact = Complex{Rational{Int}}
const physical_im = PhysicalExact(0 // 1, 1 // 1)

function physical_matrix_product(left::Matrix{PhysicalExact}, right::Matrix{PhysicalExact})
  return left * right
end

function physical_fixture()
  M0 = PhysicalExact[1 2; 0 -1]
  M1 = PhysicalExact[0 1 + physical_im; 2 -1]
  M2 = PhysicalExact[2 0; 1 physical_im]
  Mm1 = PhysicalExact[1 - physical_im 0; 1 2]
  Mm2 = PhysicalExact[-1 1; -physical_im 2]

  N0 = PhysicalExact[1 0; 2 -2]
  N1 = PhysicalExact[0 1; -1 1]
  Nm1 = PhysicalExact[2 -1; 0 -1]

  order1 = [
    physical_jump_vertex(1, :a, 0, 0, M0),
    physical_jump_vertex(1, :a, 1, 0, M1),
    physical_jump_vertex(1, :a, 2, 1, M2),
    physical_jump_vertex(1, :b, -1, 0, Mm1),
    physical_jump_vertex(1, :b, -2, -1, Mm2),
  ]
  order2 = [
    physical_drift_vertex(2, :drift, 0, 0, N0),
    physical_drift_vertex(2, :drift, 1, 0, N1),
    physical_drift_vertex(2, :drift, -1, 0, Nm1),
  ]

  zero_component = zeros(PhysicalExact, 2, 2)
  identity_component = PhysicalExact[1 0; 0 1]
  return (; vertices=[order1, order2], zero_component, identity_component)
end

@testset "physical provenance refines rather than changes the projection recurrence" begin
  fixture = physical_fixture()
  physical = compile_physical_openleg_evaluation_plan(
    fixture.vertices, 4, fixture.zero_component
  )
  plan = physical.plan
  evaluated = evaluate_openleg_evaluation_plan(
    plan;
    product=physical_matrix_product,
    inverse_weight=harmonic -> physical_im * (1 // harmonic),
    identity_component=fixture.identity_component,
  )

  collapsed_input = physical_openleg_collapse_vertices(
    fixture.vertices, fixture.zero_component
  )
  reference = openleg_bloch_reference(
    collapsed_input,
    4;
    product=physical_matrix_product,
    identity_component=fixture.identity_component,
  )

  for n in 1:4
    planned_effective = collapse_openleg_plan_component(evaluated.effective[n], plan)
    planned_wave = collapse_openleg_plan_component(evaluated.wave[n], plan)
    @test openleg_equal(planned_effective, reference.effective[n])
    @test openleg_equal(planned_wave, reference.wave[n])
  end

  same_mismatch = filter(vertex -> vertex.harmonic == 1, plan.input.vertices_by_order[1])
  @test length(same_mismatch) == 2
  @test same_mismatch[1].path != same_mismatch[2].path
  @test all(iszero(state.harmonic) for states in plan.effective_support for state in states)
  @test all(!iszero(state.harmonic) for states in plan.wave_support for state in states)

  provenance = physical.provenance
  @test any(
    vertex.system_harmonic == 1 && vertex.output_sideband == 0 && vertex.mismatch == 1 for
    vertex in provenance
  )
  @test any(
    vertex.system_harmonic == 2 && vertex.output_sideband == 1 && vertex.mismatch == 1 for
    vertex in provenance
  )
end

@testset "history path and physical Kraus-output key are distinct" begin
  first_leg = OpenLegOutputLeg(:a, 0)
  second_leg = OpenLegOutputLeg(:b, 1)
  forward = (first_leg, second_leg)
  reverse = (second_leg, first_leg)

  @test forward != reverse
  @test physical_openleg_output_key(forward) == physical_openleg_output_key(reverse)

  fixture = physical_fixture()
  physical = compile_physical_openleg_evaluation_plan(
    fixture.vertices, 5, fixture.zero_component
  )
  growth = physical_openleg_growth_counts(physical.plan)

  @test growth.raw_contributions == growth.generator_products + growth.fold_products
  @test growth.history_paths > growth.physical_output_classes
  @test growth.physical_output_coalescences > 0
  @test growth.left_right_paste_pairings >= growth.history_paths
  @test growth.recurrence_states ==
    physical.plan.counts.wave_nodes + physical.plan.counts.effective_nodes
  @test growth.dag_nodes == physical.plan.counts.residual_nodes

  # Exact order-5 structural profile for this four-label physical output support.
  @test growth.raw_contributions == 5704
  @test growth.history_paths == 1364
  @test growth.physical_output_classes == 125
  @test growth.recurrence_states == 3227
  @test growth.dag_nodes == 3227
  @test growth.generator_products == 3152
  @test growth.fold_products == 2552
  @test growth.physical_output_coalescences == 1239
  @test growth.left_right_paste_pairings == 34508
end

function physical_recycling_super(left::Matrix{PhysicalExact}, right::Matrix{PhysicalExact})
  return kron(conj.(right), left)
end

function physical_hamiltonian_super(H::Matrix{PhysicalExact})
  I2 = Matrix{PhysicalExact}(I, 2, 2)
  return -physical_im * (kron(I2, H) - kron(transpose(H), I2))
end

@testset "physical half-chain pasting reconstructs the rotating-jump RR shift" begin
  σx = PhysicalExact[0 1; 1 0]
  σy = PhysicalExact[0 -physical_im; physical_im 0]
  σz = PhysicalExact[1 0; 0 -1]
  σminus = (1 // 2) * (σx - physical_im * σy)
  σplus = (1 // 2) * (σx + physical_im * σy)
  amplitudes = Dict(0 => σz, 1 => physical_im * σminus, -1 => -physical_im * σplus)
  zero_super = zeros(PhysicalExact, 4, 4)

  rr1 = halfchain_rr_pasted_component(
    amplitudes,
    1;
    product=physical_matrix_product,
    paste=physical_recycling_super,
    minus_imaginary=(-physical_im),
    zero_component=zero_super,
  )
  rr2 = halfchain_rr_pasted_component(
    amplitudes,
    2;
    product=physical_matrix_product,
    paste=physical_recycling_super,
    minus_imaginary=(-physical_im),
    zero_component=zero_super,
  )

  @test rr1 == physical_hamiltonian_super(-σz)
  @test rr2 == physical_hamiltonian_super(-(1 // 4) * σz)
  @test rr1 + rr2 == physical_hamiltonian_super(-(5 // 4) * σz)
end

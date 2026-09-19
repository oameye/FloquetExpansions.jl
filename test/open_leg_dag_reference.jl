using Test
using LinearAlgebra: I

include(joinpath(@__DIR__, "helpers", "open_leg_reference.jl"))
include(joinpath(@__DIR__, "helpers", "bloch_evaluation_plan.jl"))
include(joinpath(@__DIR__, "helpers", "open_leg_evaluation_plan.jl"))
include(joinpath(@__DIR__, "helpers", "open_leg_hori_deprit_order3.jl"))
include(joinpath(@__DIR__, "helpers", "open_leg_bloch_normalization.jl"))

const OpenLegExactComplex = Complex{Rational{Int}}
const openleg_exact_im = OpenLegExactComplex(0 // 1, 1 // 1)

function openleg_matrix_product(
  left::Matrix{OpenLegExactComplex}, right::Matrix{OpenLegExactComplex}
)
  return left * right
end

function openleg_fixture()
  M0 = OpenLegExactComplex[1 2; 0 -1]
  M1 = OpenLegExactComplex[0 1 + openleg_exact_im; 2 -1]
  Mm1 = OpenLegExactComplex[1 - openleg_exact_im 0; 1 2]
  M2 = OpenLegExactComplex[2 0; 1 openleg_exact_im]
  Mm2 = OpenLegExactComplex[-1 1; -openleg_exact_im 2]

  N0 = OpenLegExactComplex[1 0; 2 -2]
  N1 = OpenLegExactComplex[0 1; -1 1]
  Nm1 = OpenLegExactComplex[2 -1; 0 -1]

  zero_component = zeros(OpenLegExactComplex, 2, 2)
  A1 = openleg_periodic(
    Dict((0, 1) => M0, (1, 1) => M1, (-1, 1) => Mm1, (2, 1) => M2, (-2, 1) => Mm2),
    zero_component,
  )
  A2 = openleg_periodic(Dict((0, 0) => N0, (1, 0) => N1, (-1, 0) => Nm1), zero_component)
  identity_component = OpenLegExactComplex[1 0; 0 1]
  return (; A1, A2, identity_component, zero_component)
end

function symbolic_openleg_factory(order, harmonic, grade)
  iszero(grade) && return ()
  return ntuple(
    index -> OpenLegOutputLeg(Symbol("channel_$(order)"), (harmonic, index)), grade
  )
end

@testset "open-leg DAG refines the certified Bloch recurrence" begin
  fixture = openleg_fixture()
  plan = compile_openleg_evaluation_plan(
    [fixture.A1, fixture.A2], 4; leg_factory=symbolic_openleg_factory
  )
  evaluated = evaluate_openleg_evaluation_plan(
    plan;
    product=openleg_matrix_product,
    inverse_weight=harmonic -> openleg_exact_im * (1 // harmonic),
    identity_component=fixture.identity_component,
  )
  reference = openleg_bloch_reference(
    [fixture.A1, fixture.A2],
    4;
    product=openleg_matrix_product,
    identity_component=fixture.identity_component,
  )

  for n in 1:4
    planned_effective = collapse_openleg_plan_component(evaluated.effective[n], plan)
    planned_wave = collapse_openleg_plan_component(evaluated.wave[n], plan)
    @test openleg_equal(planned_effective, reference.effective[n])
    @test openleg_equal(planned_wave, reference.wave[n])
  end

  @test any(
    state.grade == 2 && iszero(state.harmonic) for state in plan.effective_support[2]
  )
  @test all(!iszero(state.harmonic) for states in plan.wave_support for state in states)
  @test all(iszero(state.harmonic) for states in plan.effective_support for state in states)
  @test plan.counts.generator_products > 0
  @test plan.counts.fold_products > 0
  @test plan.counts.symbolic_paths > 1

  two_output_states = filter(state -> state.grade == 2, plan.effective_support[2])
  @test !isempty(two_output_states)
  @test all(length(openleg_plan_path(plan, state)) == 2 for state in two_output_states)

  primitive = [OpenLegPlanState(1, 1, 1), OpenLegPlanState(-1, 1, 1)]
  folded = [
    OpenLegPlanState(0, 1, 0), OpenLegPlanState(1, 1, 1), OpenLegPlanState(-1, 1, 1)
  ]
  @test openleg_first_return_mismatch(primitive)
  @test !openleg_first_return_mismatch(folded)
end

@testset "open-leg DAG reduces exactly to the #144 sparse evaluation schedule" begin
  fixture = openleg_fixture()
  components = Dict(harmonic => fixture.A1[harmonic, 1] for harmonic in (-2, -1, 0, 1, 2))
  grade_zero = openleg_periodic(
    Dict((harmonic, 0) => component for (harmonic, component) in components),
    fixture.zero_component,
  )

  open_plan = compile_openleg_evaluation_plan(
    [grade_zero], 5; leg_factory=(order, harmonic, grade) -> ()
  )
  open_evaluated = evaluate_openleg_evaluation_plan(
    open_plan;
    product=openleg_matrix_product,
    inverse_weight=harmonic -> openleg_exact_im * (1 // harmonic),
    identity_component=fixture.identity_component,
  )

  bloch_plan = compile_bloch_evaluation_plan(collect(keys(components)), 5)
  bloch_evaluated = evaluate_bloch_evaluation_plan(
    bloch_plan,
    components;
    product=openleg_matrix_product,
    inverse_weight=harmonic -> openleg_exact_im * (1 // harmonic),
    zero_component=fixture.zero_component,
  )

  for n in 1:5
    collapsed = collapse_openleg_plan_component(open_evaluated.effective[n], open_plan)
    @test collapsed[0, 0] == bloch_evaluated.effective[n]
  end

  for n in 1:4
    collapsed = collapse_openleg_plan_component(open_evaluated.wave[n], open_plan)
    for harmonic in keys(bloch_evaluated.wave[n])
      @test collapsed[harmonic, 0] == bloch_evaluated.wave[n][harmonic]
    end
  end
end

@testset "Hori-Deprit and Bloch-Feshbach agree through the one-output cubic sector" begin
  fixture = openleg_fixture()
  bloch = openleg_bloch_reference(
    [fixture.A1, fixture.A2],
    3;
    product=openleg_matrix_product,
    identity_component=fixture.identity_component,
  )
  hori = openleg_hori_deprit_order3(fixture.A1, fixture.A2; product=openleg_matrix_product)

  @test bloch.effective[2][0, 0] == hori.effective2[0, 0]
  @test bloch.effective[2][0, 2] == hori.effective2[0, 2]
  @test bloch.effective[3][0, 1] == hori.effective3[0, 1]

  # Raw Bloch and zero-average Lie gauges first differ in the three-output cubic sector.
  @test bloch.effective[3][0, 3] != hori.effective3[0, 3]

  static2 = openleg_bloch_static_factor2(bloch; product=openleg_matrix_product)
  @test openleg_grades(static2) == [2]
  canonical_b3 = openleg_bloch_canonical_effective3(bloch; product=openleg_matrix_product)
  @test openleg_equal(canonical_b3, hori.effective3)
end

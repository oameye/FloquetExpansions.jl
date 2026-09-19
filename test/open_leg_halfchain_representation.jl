using Test

include(joinpath(@__DIR__, "helpers", "open_leg_reference.jl"))
include(joinpath(@__DIR__, "helpers", "open_leg_evaluation_plan.jl"))
include(joinpath(@__DIR__, "helpers", "stinespring_normalization_reference.jl"))
include(joinpath(@__DIR__, "helpers", "open_leg_physical_provenance.jl"))
include(joinpath(@__DIR__, "helpers", "open_leg_halfchain_representation.jl"))

const HalfChainExact = Complex{Rational{Int}}
const halfchain_im = HalfChainExact(0 // 1, 1 // 1)

function halfchain_matrix_product(
  left::Matrix{HalfChainExact}, right::Matrix{HalfChainExact}
)
  return left * right
end

function halfchain_physical_fixture()
  M0 = HalfChainExact[1 2; 0 -1]
  M1 = HalfChainExact[0 1 + halfchain_im; 2 -1]
  M2 = HalfChainExact[2 0; 1 halfchain_im]
  Mm1 = HalfChainExact[1 - halfchain_im 0; 1 2]
  Mm2 = HalfChainExact[-1 1; -halfchain_im 2]

  N0 = HalfChainExact[1 0; 2 -2]
  N1 = HalfChainExact[0 1; -1 1]
  Nm1 = HalfChainExact[2 -1; 0 -1]

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

  zero_component = zeros(HalfChainExact, 2, 2)
  identity_component = HalfChainExact[1 0; 0 1]
  return (; vertices=[order1, order2], zero_component, identity_component)
end

function sparse_scalar_fixture(output_labels::Int)
  output_labels in (2, 3) || throw(ArgumentError("fixture supports two or three labels"))
  order1 = [
    physical_jump_vertex(1, :a, 0, 0, 1 // 1),
    physical_jump_vertex(1, :a, 1, 0, 2 // 1),
    physical_jump_vertex(1, :b, -1, 0, 3 // 1),
  ]
  if output_labels == 3
    push!(order1, physical_jump_vertex(1, :b, -2, -1, 4 // 1))
  end
  order2 = [
    physical_drift_vertex(2, :drift, 0, 0, 1 // 1),
    physical_drift_vertex(2, :drift, 1, 0, 2 // 1),
    physical_drift_vertex(2, :drift, -1, 0, 3 // 1),
  ]
  return (; vertices=[order1, order2], zero_component=0 // 1, identity_component=1 // 1)
end

@testset "persistent paths preserve exact ordered-history identity" begin
  fixture = halfchain_physical_fixture()
  flat =
    compile_physical_openleg_evaluation_plan(fixture.vertices, 5, fixture.zero_component).plan
  representation = physical_openleg_representation_counts(flat)

  @test representation.flat_paths == 1364
  @test representation.flat_leg_cells == 6372
  @test representation.persistent_nodes == 1364
  @test representation.persistent_leg_atoms == 4
  @test representation.physical_output_classes == 125
  @test representation.physical_output_key_cells == 504

  @test representation.persistent_nodes < representation.flat_leg_cells
  @test representation.physical_output_classes < representation.flat_paths
  @test representation.physical_output_key_cells < representation.flat_leg_cells
end

@testset "physical output quotient respects composition" begin
  fixture = halfchain_physical_fixture()
  flat =
    compile_physical_openleg_evaluation_plan(fixture.vertices, 3, fixture.zero_component).plan
  outputs = PhysicalOpenLegOutputRegistry()
  output_ids = [physical_openleg_intern_output!(outputs, path) for path in flat.input.paths]

  for (left, left_id) in zip(flat.input.paths, output_ids),
    (right, right_id) in zip(flat.input.paths, output_ids)

    composed = physical_openleg_compose_output!(outputs, left_id, right_id)
    @test outputs.keys[composed] == physical_openleg_output_key((left..., right...))
  end

  for states in (flat.wave_support..., flat.effective_support...), state in states
    @test state.grade == length(physical_openleg_output_key(flat.input.paths[state.path]))
  end
  @test fieldnames(PhysicalOutputState) == (:harmonic, :output)
end

@testset "physical output sector is a sufficient recurrence state" begin
  fixture = halfchain_physical_fixture()
  flat =
    compile_physical_openleg_evaluation_plan(fixture.vertices, 4, fixture.zero_component).plan
  flat_evaluated = evaluate_openleg_evaluation_plan(
    flat;
    product=halfchain_matrix_product,
    inverse_weight=harmonic -> halfchain_im * (1 // harmonic),
    identity_component=fixture.identity_component,
  )
  compact = evaluate_physical_output_recurrence(
    fixture.vertices,
    4,
    fixture.zero_component;
    product=halfchain_matrix_product,
    inverse_weight=harmonic -> halfchain_im * (1 // harmonic),
    identity_component=fixture.identity_component,
  )

  for n in 1:4
    @test materialize_compact_physical_outputs(
      compact.wave[n], compact.outputs, fixture.zero_component
    ) == coalesce_flat_physical_outputs(flat_evaluated.wave[n], flat)
    @test materialize_compact_physical_outputs(
      compact.effective[n], compact.outputs, fixture.zero_component
    ) == coalesce_flat_physical_outputs(flat_evaluated.effective[n], flat)
  end

  @test all(!iszero(state.harmonic) for states in compact.wave for state in keys(states))
  @test all(
    iszero(state.harmonic) for states in compact.effective for state in keys(states)
  )
end

@testset "sufficient-state recurrence removes ordered-history state growth" begin
  fixture = halfchain_physical_fixture()
  flat =
    compile_physical_openleg_evaluation_plan(fixture.vertices, 5, fixture.zero_component).plan
  flat_growth = physical_openleg_growth_counts(flat)
  compact = evaluate_physical_output_recurrence(
    fixture.vertices,
    5,
    fixture.zero_component;
    product=halfchain_matrix_product,
    inverse_weight=harmonic -> halfchain_im * (1 // harmonic),
    identity_component=fixture.identity_component,
  )

  @test compact.counts.residual_states == 399
  @test compact.counts.wave_states == 334
  @test compact.counts.effective_states == 65
  @test compact.counts.generator_products == 893
  @test compact.counts.fold_products == 627
  @test compact.counts.physical_output_classes == 125

  @test compact.counts.residual_states < flat_growth.recurrence_states
  @test compact.counts.generator_products + compact.counts.fold_products <
    flat_growth.raw_contributions
end

@testset "compression persists across sparse physical supports" begin
  for fixture in (sparse_scalar_fixture(2), sparse_scalar_fixture(3))
    flat =
      compile_physical_openleg_evaluation_plan(fixture.vertices, 5, fixture.zero_component).plan
    flat_growth = physical_openleg_growth_counts(flat)
    representation = physical_openleg_representation_counts(flat)
    compact = evaluate_physical_output_recurrence(
      fixture.vertices,
      5,
      fixture.zero_component;
      product=(*),
      inverse_weight=harmonic -> 1 // harmonic,
      identity_component=fixture.identity_component,
    )

    @test representation.persistent_nodes <= representation.flat_paths
    @test representation.persistent_nodes < representation.flat_leg_cells
    @test representation.physical_output_classes < representation.flat_paths
    @test compact.counts.residual_states < flat_growth.recurrence_states
    @test compact.counts.generator_products + compact.counts.fold_products <
      flat_growth.raw_contributions
  end
end

using Test
using FloquetExpansions

include(joinpath(@__DIR__, "helpers", "ck_graded_output_algebra.jl"))

const CKExact = Complex{Rational{Int}}
const ck_im = CKExact(0 // 1, 1 // 1)

function ck_physical_fixture()
  M0 = CKExact[1 2; 0 -1]
  M1 = CKExact[0 1 + ck_im; 2 -1]
  M2 = CKExact[2 0; 1 ck_im]
  Mm1 = CKExact[1 - ck_im 0; 1 2]
  Mm2 = CKExact[-1 1; -ck_im 2]

  N0 = CKExact[1 0; 2 -2]
  N1 = CKExact[0 1; -1 1]
  Nm1 = CKExact[2 -1; 0 -1]

  order1 = [
    ck_jump_vertex(1, :a, 0, 0, M0),
    ck_jump_vertex(1, :a, 1, 0, M1),
    ck_jump_vertex(1, :a, 2, 1, M2),
    ck_jump_vertex(1, :b, -1, 0, Mm1),
    ck_jump_vertex(1, :b, -2, -1, Mm2),
  ]
  order2 = [
    ck_drift_vertex(2, :drift, 0, 0, N0),
    ck_drift_vertex(2, :drift, 1, 0, N1),
    ck_drift_vertex(2, :drift, -1, 0, Nm1),
  ]

  return (
    vertices=[order1, order2],
    zero_component=zeros(CKExact, 2, 2),
    identity_component=CKExact[1 0; 0 1],
  )
end

@testset "production graded output algebra is concrete and internal" begin
  left = @inferred FloquetExpansions.graded_output_series(Int, 0 // 1, 3)
  right = @inferred FloquetExpansions.graded_output_series(Int, 0 // 1, 3)
  @test @inferred(FloquetExpansions.graded_output_accumulate!(left, 1, 1, 2 // 1)) === left
  @test @inferred(FloquetExpansions.graded_output_accumulate!(right, 1, 2, 3 // 1)) ===
    right

  compose_output =
    (left_key, right_key) -> FloquetExpansions.OutputSectorComposition(left_key + right_key)
  algebra = FloquetExpansions.GradedOutputProduct(compose_output, *)
  product = @inferred algebra(left, right)
  @test product isa FloquetExpansions.GradedOutputSeries{Int,Rational{Int}}
  @test product.coefficients[3] == Dict(3 => 6 // 1)
  @test @inferred((2 // 1) * left).coefficients[2] == Dict(1 => 4 // 1)

  weighted_compose =
    (left_key, right_key) ->
      FloquetExpansions.OutputSectorComposition(left_key + right_key, -1)
  weighted_product = @inferred FloquetExpansions.GradedOutputProduct(weighted_compose, *)(
    left, right
  )
  @test weighted_product.coefficients[3] == Dict(3 => -6 // 1)

  @test_throws ArgumentError left + FloquetExpansions.graded_output_series(Int, 0 // 1, 2)
  @test !isdefined(Main, :GradedOutputSeries)
end

@testset "graded CK amplitudes reuse the shared Bloch projection core" begin
  fixture = ck_physical_fixture()
  cutoff = 5
  components = ck_series_components(fixture.vertices, fixture.zero_component, cutoff)
  plan = FloquetExpansions.compile_bloch_projection_plan(keys(components), cutoff + 1)
  zero_series = ck_zero_series(fixture.zero_component, cutoff)
  series_products = Ref(0)
  series_product = ck_counted_series_product(series_products)

  bloch = FloquetExpansions.evaluate_bloch_projection_plan(
    plan,
    components;
    product=series_product,
    inverse_weight=harmonic -> ck_im * (1 // harmonic),
    zero_component=zero_series,
  )
  direct = ck_direct_tiered_recurrence(
    fixture.vertices,
    cutoff,
    fixture.zero_component,
    fixture.identity_component;
    inverse_weight=harmonic -> ck_im * (1 // harmonic),
  )

  for degree in 1:cutoff
    @test ck_materialize_bloch_wave(bloch, degree, fixture.zero_component) ==
      direct.wave[degree]
    @test ck_materialize_bloch_effective(bloch, degree, fixture.zero_component) ==
      direct.effective[degree]
  end

  @test direct.products == 1520
  @test series_products[] == 1512
  @test series_products[] == direct.products - sum(length, fixture.vertices)
end

@testset "single-tier specialization reduces to ordinary Bloch evaluation" begin
  fixture = ck_physical_fixture()
  zero_component = fixture.zero_component
  components = Dict(
    -1 => fixture.vertices[1][4].value,
    0 => fixture.vertices[1][1].value,
    1 => fixture.vertices[1][2].value,
  )
  cutoff = 5
  plan = FloquetExpansions.compile_bloch_projection_plan(keys(components), cutoff)

  plain = FloquetExpansions.evaluate_bloch_projection_plan(
    plan,
    components;
    product=(*),
    inverse_weight=harmonic -> ck_im * (1 // harmonic),
    zero_component,
  )

  series_components = ck_trivial_series_components(components, zero_component, cutoff)
  series_products = Ref(0)
  series_product = ck_counted_series_product(series_products)
  series = FloquetExpansions.evaluate_bloch_projection_plan(
    plan,
    series_components;
    product=series_product,
    inverse_weight=harmonic -> ck_im * (1 // harmonic),
    zero_component=ck_zero_series(zero_component, cutoff),
  )

  for degree in 1:cutoff
    effective = ck_materialize_bloch_effective(series, degree, zero_component)
    expected = plain.effective[degree]
    expected_state = CKDirectState(0, ())
    if expected == zero_component
      @test isempty(effective)
    else
      @test effective == Dict(expected_state => expected)
    end
  end

  for degree in 1:(cutoff - 1)
    wave = ck_materialize_bloch_wave(series, degree, zero_component)
    expected = Dict(
      CKDirectState(harmonic, ()) => value for (harmonic, value) in plain.wave[degree]
    )
    @test wave == expected
  end
end

struct CKTwoHarmonic
  first::Int
  second::Int
end

Base.zero(::CKTwoHarmonic) = CKTwoHarmonic(0, 0)
Base.iszero(harmonic::CKTwoHarmonic) = iszero(harmonic.first) && iszero(harmonic.second)
function Base.:+(left::CKTwoHarmonic, right::CKTwoHarmonic)
  return CKTwoHarmonic(left.first + right.first, left.second + right.second)
end

@testset "graded coefficient algebra preserves generic harmonic labels" begin
  zero_harmonic = CKTwoHarmonic(0, 0)
  minus = CKTwoHarmonic(-1, 0)
  plus = CKTwoHarmonic(1, 0)
  components = Dict(minus => 2 // 1, zero_harmonic => 3 // 1, plus => 5 // 1)
  cutoff = 4
  plan = FloquetExpansions.compile_bloch_projection_plan(keys(components), cutoff)
  plain = FloquetExpansions.evaluate_bloch_projection_plan(
    plan,
    components;
    product=(*),
    inverse_weight=harmonic -> 1 // harmonic.first,
    zero_component=0 // 1,
  )

  series_components = ck_trivial_series_components(components, 0 // 1, cutoff)
  series_products = Ref(0)
  series_product = ck_counted_series_product(series_products)
  series = FloquetExpansions.evaluate_bloch_projection_plan(
    plan,
    series_components;
    product=series_product,
    inverse_weight=harmonic -> 1 // harmonic.first,
    zero_component=ck_zero_series(0 // 1, cutoff),
  )

  for degree in 1:cutoff
    coefficient = get(series.effective[degree].coefficients[degree + 1], (), 0 // 1)
    @test coefficient == plain.effective[degree]
  end
  for degree in 1:(cutoff - 1), harmonic in keys(plain.wave[degree])
    coefficient = get(series.wave[degree][harmonic].coefficients[degree + 1], (), 0 // 1)
    @test coefficient == plain.wave[degree][harmonic]
  end

  @test plan isa FloquetExpansions.BlochProjectionPlan{CKTwoHarmonic}
end

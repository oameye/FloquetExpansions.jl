using Test
using FloquetExpansions
using LinearAlgebra: I

struct BlochOrderTestSeries{T}
  components::Dict{Int,T}
  zero_component::T
end

function bloch_order_test_series(
  components::AbstractDict{Int,T}, zero_component::T
) where {T}
  return BlochOrderTestSeries(
    Dict(harmonic => value for (harmonic, value) in components if value != zero_component),
    zero_component,
  )
end

function Base.:(==)(left::BlochOrderTestSeries, right::BlochOrderTestSeries)
  return left.zero_component == right.zero_component && left.components == right.components
end

function Base.:+(left::BlochOrderTestSeries{T}, right::BlochOrderTestSeries{T}) where {T}
  result = copy(left.components)
  for (harmonic, value) in right.components
    updated = get(result, harmonic, left.zero_component) + value
    if updated == left.zero_component
      haskey(result, harmonic) && delete!(result, harmonic)
    else
      result[harmonic] = updated
    end
  end
  return BlochOrderTestSeries(result, left.zero_component)
end

function Base.:-(left::BlochOrderTestSeries{T}, right::BlochOrderTestSeries{T}) where {T}
  result = copy(left.components)
  for (harmonic, value) in right.components
    updated = get(result, harmonic, left.zero_component) - value
    if updated == left.zero_component
      haskey(result, harmonic) && delete!(result, harmonic)
    else
      result[harmonic] = updated
    end
  end
  return BlochOrderTestSeries(result, left.zero_component)
end

function bloch_order_test_product(
  left::BlochOrderTestSeries{T}, right::BlochOrderTestSeries{T}
) where {T}
  result = Dict{Int,T}()
  for (left_harmonic, left_value) in left.components,
    (right_harmonic, right_value) in right.components

    harmonic = left_harmonic + right_harmonic
    updated = get(result, harmonic, left.zero_component) + left_value * right_value
    if updated == left.zero_component
      haskey(result, harmonic) && delete!(result, harmonic)
    else
      result[harmonic] = updated
    end
  end
  return BlochOrderTestSeries(result, left.zero_component)
end

function bloch_order_test_project_model(series::BlochOrderTestSeries{T}) where {T}
  value = get(series.components, 0, series.zero_component)
  return if value == series.zero_component
    BlochOrderTestSeries(Dict{Int,T}(), series.zero_component)
  else
    BlochOrderTestSeries(Dict(0 => value), series.zero_component)
  end
end

function bloch_order_test_solve_complement(series::BlochOrderTestSeries{T}) where {T}
  result = Dict{Int,T}()
  for (harmonic, value) in series.components
    iszero(harmonic) && continue
    result[harmonic] = (1 // harmonic) * value
  end
  return BlochOrderTestSeries(result, series.zero_component)
end

function bloch_order_test_value(series::BlochOrderTestSeries, harmonic::Int)
  return get(series.components, harmonic, series.zero_component)
end

function bloch_order_test_operations()
  return FloquetExpansions.BlochOrderOperations(
    bloch_order_test_product,
    bloch_order_test_project_model,
    bloch_order_test_solve_complement,
  )
end

@testset "generic Bloch order recurrence is concrete and internal" begin
  zero_component = zeros(Rational{Int}, 2, 2)
  identity_component = Matrix{Rational{Int}}(I, 2, 2)
  A1 = bloch_order_test_series(
    Dict(
      -1 => Rational{Int}[1 2; 0 -1],
      0 => Rational{Int}[0 1; 1 0],
      1 => Rational{Int}[2 -1; 1 1],
    ),
    zero_component,
  )
  identity_series = BlochOrderTestSeries(Dict(0 => identity_component), zero_component)
  zero_series = BlochOrderTestSeries(Dict{Int,Matrix{Rational{Int}}}(), zero_component)
  operations = bloch_order_test_operations()

  result = @inferred FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1], 5, identity_series, zero_series, operations
  )

  @test result isa FloquetExpansions.BlochOrderRecurrenceResult{typeof(A1)}
  @test length(result.effective) == 5
  @test length(result.wave) == 4
  @test result.products == 15
  @test_throws ArgumentError FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1], 0, identity_series, zero_series, operations
  )
  @test_throws ArgumentError FloquetExpansions.evaluate_bloch_order_recurrence(
    typeof(A1)[], 1, identity_series, zero_series, operations
  )
  @test !isdefined(Main, :evaluate_bloch_order_recurrence)
  @test !isdefined(Main, :BlochOrderOperations)
end

@testset "generic recurrence exactly reproduces the finite Fourier projection core" begin
  zero_component = zeros(Rational{Int}, 2, 2)
  identity_component = Matrix{Rational{Int}}(I, 2, 2)
  components = Dict(
    -2 => Rational{Int}[1 0; 2 -1],
    -1 => Rational{Int}[0 1; -1 2],
    0 => Rational{Int}[1 2; 3 -2],
    1 => Rational{Int}[2 -1; 1 0],
    2 => Rational{Int}[-1 1; 2 1],
  )
  A1 = bloch_order_test_series(components, zero_component)
  identity_series = BlochOrderTestSeries(Dict(0 => identity_component), zero_component)
  zero_series = BlochOrderTestSeries(Dict{Int,Matrix{Rational{Int}}}(), zero_component)
  operations = bloch_order_test_operations()
  order = 5

  generic = FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1], order, identity_series, zero_series, operations
  )

  plan = FloquetExpansions.compile_bloch_projection_plan(keys(components), order)
  planned = FloquetExpansions.evaluate_bloch_projection_plan(
    plan, components; product=(*), inverse_weight=harmonic -> 1 // harmonic, zero_component
  )

  for n in 1:order
    @test bloch_order_test_value(generic.effective[n], 0) == planned.effective[n]
  end
  for n in 1:(order - 1)
    harmonics = union(keys(generic.wave[n].components), keys(planned.wave[n]))
    for harmonic in harmonics
      @test bloch_order_test_value(generic.wave[n], harmonic) ==
        get(planned.wave[n], harmonic, zero_component)
    end
  end
end

@testset "tiered arbitrary-order recurrence is prefix consistent" begin
  zero_component = zeros(Rational{Int}, 2, 2)
  identity_component = Matrix{Rational{Int}}(I, 2, 2)
  A1 = bloch_order_test_series(
    Dict(
      -1 => Rational{Int}[1 1; 0 -1],
      0 => Rational{Int}[0 1; 1 1],
      1 => Rational{Int}[2 -1; 1 0],
    ),
    zero_component,
  )
  A2 = bloch_order_test_series(
    Dict(0 => Rational{Int}[1 -1; 2 0], 1 => Rational{Int}[0 2; -1 1]), zero_component
  )
  identity_series = BlochOrderTestSeries(Dict(0 => identity_component), zero_component)
  zero_series = BlochOrderTestSeries(Dict{Int,Matrix{Rational{Int}}}(), zero_component)
  operations = bloch_order_test_operations()

  order4 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1, A2], 4, identity_series, zero_series, operations
  )
  order5 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1, A2], 5, identity_series, zero_series, operations
  )
  one_tier5 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [A1], 5, identity_series, zero_series, operations
  )

  @test order4.effective == order5.effective[1:4]
  @test order4.wave == order5.wave[1:3]
  @test one_tier5.products == 15
  @test order4.products == 13
  @test order5.products == 19
end

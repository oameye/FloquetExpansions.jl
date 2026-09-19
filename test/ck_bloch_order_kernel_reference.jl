using Test
using LinearAlgebra: I

include(joinpath(@__DIR__, "..", "src", "bloch_projection.jl"))
include(joinpath(@__DIR__, "helpers", "bloch_order_kernel_reference.jl"))

const CKOrderExact = Complex{Rational{Int}}
const ck_order_im = CKOrderExact(0 // 1, 1 // 1)

function ck_order_fixture()
  zero_component = zeros(CKOrderExact, 2, 2)
  identity_component = Matrix{CKOrderExact}(I, 2, 2)
  first = Dict(
    -1 => CKOrderExact[0 1; 2 -1],
    0 => CKOrderExact[1 2; -1 0],
    1 => CKOrderExact[2 -1; 1 1],
  )
  second = Dict(
    -1 => CKOrderExact[1 0; -1 2], 0 => CKOrderExact[0 1; 1 -1], 1 => CKOrderExact[2 1; 0 1]
  )
  return (; zero_component, identity_component, first, second)
end

@testset "order-only Bloch recurrence reproduces #148 exactly on finite Fourier data" begin
  fixture = ck_order_fixture()
  order = 5
  inverse_weight(harmonic) = ck_order_im * (1 // harmonic)

  plan = compile_bloch_projection_plan(keys(fixture.first), order)
  discrete = evaluate_bloch_projection_plan(
    plan, fixture.first; product=(*), inverse_weight, zero_component=fixture.zero_component
  )

  generator = kernel_fourier_series(fixture.first, fixture.zero_component)
  identity_series = kernel_fourier_series(
    Dict(0 => fixture.identity_component), fixture.zero_component
  )
  zero_series = kernel_fourier_series(
    Dict{Int,Matrix{CKOrderExact}}(), fixture.zero_component
  )
  packed = bloch_order_recurrence(
    [generator],
    order,
    identity_series,
    zero_series;
    product=kernel_fourier_product,
    project_model=kernel_fourier_project_model,
    solve_complement=series -> kernel_fourier_solve_complement(series, inverse_weight),
  )

  for n in 1:order
    @test get(packed.effective[n].components, 0, fixture.zero_component) ==
      discrete.effective[n]
    @test all(iszero, keys(packed.effective[n].components))
  end
  for n in 1:(order - 1)
    @test packed.wave[n].components == discrete.wave[n]
  end

  @test packed.products == order * (order + 1) ÷ 2
  @test !isempty(kernel_fourier_project_model(generator).components)
  @test !isempty(kernel_fourier_solve_complement(generator, inverse_weight).components)
end

@testset "tiered kernel recurrence is prefix-consistent without harmonic-key routing" begin
  fixture = ck_order_fixture()
  inverse_weight(harmonic) = ck_order_im * (1 // harmonic)
  first = kernel_fourier_series(fixture.first, fixture.zero_component)
  second = kernel_fourier_series(fixture.second, fixture.zero_component)
  identity_series = kernel_fourier_series(
    Dict(0 => fixture.identity_component), fixture.zero_component
  )
  zero_series = kernel_fourier_series(
    Dict{Int,Matrix{CKOrderExact}}(), fixture.zero_component
  )

  solve(series) = kernel_fourier_solve_complement(series, inverse_weight)
  short = bloch_order_recurrence(
    [first, second],
    4,
    identity_series,
    zero_series;
    product=kernel_fourier_product,
    project_model=kernel_fourier_project_model,
    solve_complement=solve,
  )
  long = bloch_order_recurrence(
    [first, second],
    5,
    identity_series,
    zero_series;
    product=kernel_fourier_product,
    project_model=kernel_fourier_project_model,
    solve_complement=solve,
  )

  for n in 1:4
    @test short.effective[n].components == long.effective[n].components
  end
  for n in 1:3
    @test short.wave[n].components == long.wave[n].components
  end

  @test short.products == 13
  @test long.products == 19
end

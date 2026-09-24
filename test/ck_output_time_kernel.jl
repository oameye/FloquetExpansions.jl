using Test
using FloquetExpansions

const CKOutputTimeExact = Complex{Rational{Int}}
const ck_output_time_im = CKOutputTimeExact(0 // 1, 1 // 1)

function ck_output_time_terms(polynomial)
  return polynomial.terms
end

@testset "arbitrary-order homological kernels reproduce periodic Bernoulli polynomials" begin
  one_rational = 1 // 1
  expected = [
    Dict((0, 1) => 1 // 1, (1, 0) => -1 // 2),
    Dict((0, 2) => -1 // 2, (1, 1) => 1 // 2, (2, 0) => -1 // 12),
    Dict((0, 3) => 1 // 6, (1, 2) => -1 // 4, (2, 1) => 1 // 12),
    Dict((0, 4) => -1 // 24, (1, 3) => 1 // 12, (2, 2) => -1 // 24, (4, 0) => 1 // 720),
  ]

  kernels = [FloquetExpansions.ck_homological_kernel(order, one_rational) for order in 1:4]
  for order in 1:4
    @test ck_output_time_terms(kernels[order]) == expected[order]
    @test isempty(FloquetExpansions.ck_homological_average(kernels[order]))
  end

  for order in 1:3
    derivative = FloquetExpansions.ck_homological_derivative(kernels[order + 1])
    negative = FloquetExpansions.ck_homological_scale(kernels[order], -one_rational)
    @test derivative.terms == negative.terms
  end
end

@testset "distinct Q-resolvent poles require finite phase corrections" begin
  poles = [2, 5]
  kernel = FloquetExpansions.ck_resolvent_time_kernel(poles, ck_output_time_im)

  @test kernel.homological == [
    FloquetExpansions.CKHomologicalComponent(2, 1, ck_output_time_im / 3),
    FloquetExpansions.CKHomologicalComponent(5, 1, -ck_output_time_im / 3),
  ]
  @test kernel.phases == [
    FloquetExpansions.CKPhaseComponent(2, CKOutputTimeExact(-1 // 9)),
    FloquetExpansions.CKPhaseComponent(5, CKOutputTimeExact(-1 // 9)),
  ]

  for sideband in -8:8
    from_time = FloquetExpansions.ck_resolvent_time_fourier_coefficient(
      kernel, sideband, ck_output_time_im
    )
    direct = FloquetExpansions.ck_resolvent_product_fourier_coefficient(
      poles, sideband, ck_output_time_im
    )
    @test from_time == direct
  end
end

@testset "repeated and distinct poles reduce exactly at arbitrary multiplicity" begin
  poles = [2, 2, 5]
  kernel = FloquetExpansions.ck_resolvent_time_kernel(poles, ck_output_time_im)

  @test kernel.homological == [
    FloquetExpansions.CKHomologicalComponent(2, 1, CKOutputTimeExact(1 // 9)),
    FloquetExpansions.CKHomologicalComponent(2, 2, ck_output_time_im / 3),
    FloquetExpansions.CKHomologicalComponent(5, 1, CKOutputTimeExact(-1 // 9)),
  ]
  @test kernel.phases == [
    FloquetExpansions.CKPhaseComponent(2, ck_output_time_im / 27),
    FloquetExpansions.CKPhaseComponent(5, 2 * ck_output_time_im / 27),
  ]

  for sideband in -8:8
    @test FloquetExpansions.ck_resolvent_time_fourier_coefficient(
      kernel, sideband, ck_output_time_im
    ) == FloquetExpansions.ck_resolvent_product_fourier_coefficient(
      poles, sideband, ck_output_time_im
    )
  end

  repeated = [3, 3, 3, 3]
  repeated_kernel = FloquetExpansions.ck_resolvent_time_kernel(repeated, ck_output_time_im)
  @test repeated_kernel.homological ==
    [FloquetExpansions.CKHomologicalComponent(3, 4, one(CKOutputTimeExact))]
  @test isempty(repeated_kernel.phases)
  for sideband in -8:8
    @test FloquetExpansions.ck_resolvent_time_fourier_coefficient(
      repeated_kernel, sideband, ck_output_time_im
    ) == FloquetExpansions.ck_resolvent_product_fourier_coefficient(
      repeated, sideband, ck_output_time_im
    )
  end
end

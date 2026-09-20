using Test
using FloquetExpansions

const CKStaticTriangleExact = Complex{Rational{Int}}
const ck_static_triangle_im = CKStaticTriangleExact(0 // 1, 1 // 1)

function ck_static_triangle_key(output_number::Int)
  output_number >= 0 || throw(ArgumentError("output number must be nonnegative"))
  iszero(output_number) &&
    return FloquetExpansions.CKOutputKernelKey(
      Int[],
      FloquetExpansions.CKOutputBlock[],
      FloquetExpansions.CKOutputConstraint{Int}[],
      FloquetExpansions.CKOutputConstraint{Int}[],
    )

  return FloquetExpansions.CKOutputKernelKey(
    fill(1, output_number),
    [FloquetExpansions.CKOutputBlock(0, output_number, 0, output_number)],
    FloquetExpansions.CKOutputConstraint{Int}[
      FloquetExpansions.CKOutputConstraint(1, output_stop, 0) for
      output_stop in 1:output_number
    ],
    FloquetExpansions.CKOutputConstraint{Int}[],
  )
end

function ck_static_triangle_amplitudes(depth::Int)
  depth >= 1 || throw(ArgumentError("triangle depth must be positive"))
  max_order = 2 * depth
  zero_component = zero(CKStaticTriangleExact)
  amplitudes = FloquetExpansions.CKOutputPeriodPolynomial{Int,CKStaticTriangleExact}[]

  for perturbative_order in 0:max_order
    coefficients = FloquetExpansions.CKOutputKernel{Int,CKStaticTriangleExact}[
      FloquetExpansions.CKOutputKernel(
        Dict{FloquetExpansions.CKOutputKernelKey{Int},CKStaticTriangleExact}(),
        zero_component,
      ) for _ in 0:depth
    ]

    for output_number in 0:depth
      for drift_order in 0:(depth - output_number)
        output_number + 2 * drift_order == perturbative_order || continue
        period_power = output_number + drift_order
        value = CKStaticTriangleExact((-1 // 2)^drift_order / factorial(drift_order))
        key = ck_static_triangle_key(output_number)
        coefficients[period_power + 1].terms[key] = value
      end
    end
    push!(amplitudes, FloquetExpansions.CKOutputPeriodPolynomial(coefficients))
  end
  return amplitudes
end

ck_static_triangle_metric_pair(left, right) = conj(left) * right

function ck_static_triangle_amplitude_equal(left, right)
  length(left.coefficients) == length(right.coefficients) || return false
  return all(
    left_coefficient.zero_component == right_coefficient.zero_component &&
    left_coefficient.terms == right_coefficient.terms for
    (left_coefficient, right_coefficient) in zip(left.coefficients, right.coefficients)
  )
end

@testset "production CK static triangle has the #201 TP defect order" begin
  for generator_order in 0:2
    depth = generator_order + 1
    amplitudes = ck_static_triangle_amplitudes(depth)
    defect_order = 2 * depth + 2
    metric = FloquetExpansions.ck_output_metric_series(
      amplitudes,
      defect_order,
      ck_static_triangle_im,
      zero(CKStaticTriangleExact),
      ck_static_triangle_metric_pair,
    )

    @test FloquetExpansions.ck_output_pairing_period_terms(metric.coefficients[1]) ==
      Dict(0 => one(CKStaticTriangleExact))
    @test all(isempty(metric.coefficients[order + 1].terms) for order in 1:(defect_order - 1))
    @test !isempty(metric.coefficients[defect_order + 1].terms)
    @test FloquetExpansions.ck_output_pairing_phase_support(
      metric.coefficients[defect_order + 1]
    ) == [0]

    normalization = FloquetExpansions.ck_output_metric_inverse_sqrt_series(
      metric, defect_order, one(CKStaticTriangleExact)
    )
    @test all(
      isempty(normalization.coefficients[order + 1].terms) for
      order in 1:(defect_order - 1)
    )
    @test !isempty(normalization.coefficients[defect_order + 1].terms)

    normalized_amplitudes = FloquetExpansions.ck_output_right_normalize_series(
      amplitudes, normalization, defect_order, zero(CKStaticTriangleExact)
    )
    @test all(
      ck_static_triangle_amplitude_equal(normalized_amplitudes[index], amplitudes[index]) for
      index in eachindex(amplitudes)
    )

    normalized_metric = FloquetExpansions.ck_output_metric_series(
      normalized_amplitudes,
      defect_order,
      ck_static_triangle_im,
      zero(CKStaticTriangleExact),
      ck_static_triangle_metric_pair,
    )
    @test FloquetExpansions.ck_output_pairing_period_terms(
      normalized_metric.coefficients[1]
    ) == Dict(0 => one(CKStaticTriangleExact))
    @test all(
      isempty(normalized_metric.coefficients[order + 1].terms) for order in 1:defect_order
    )
  end
end

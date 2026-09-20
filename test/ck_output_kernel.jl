using Test
using FloquetExpansions

const CKOutputExact = Complex{Rational{Int}}
const ck_output_im = CKOutputExact(0 // 1, 1 // 1)

function ck_output_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_output_im / mismatch
end

@testset "endpoint constraints reduce to finite cumulative-output spans" begin
  zero_component = zero(CKOutputExact)
  value = CKOutputExact(3 // 2, -2 // 3)
  vertices = [
    FloquetExpansions.ck_jump_vertex(1, 2),
    FloquetExpansions.ck_drift_vertex(3),
    FloquetExpansions.ck_jump_vertex(2, -1),
  ]
  endpoint_key = FloquetExpansions.CKEndpointKey(
    vertices,
    [FloquetExpansions.CKEndpointConstraint(0, 3)],
    [
      FloquetExpansions.CKEndpointConstraint(0, 1),
      FloquetExpansions.CKEndpointConstraint(0, 2),
    ],
  )
  endpoint = FloquetExpansions.CKEndpointKernel(Dict(endpoint_key => value), zero_component)
  output = FloquetExpansions.ck_output_kernel(endpoint)

  @test length(output.terms) == 1
  output_key = only(keys(output.terms))
  @test output_key.output_channels == [1, 2]
  @test output_key.blocks == [FloquetExpansions.CKOutputBlock(0, 3, 0, 2)]
  @test output_key.model_constraints == [FloquetExpansions.CKOutputConstraint(1, 2, 4)]
  @test output_key.resolvent_constraints == [
    FloquetExpansions.CKOutputConstraint(1, 1, 2),
    FloquetExpansions.CKOutputConstraint(1, 1, 5),
  ]
  @test FloquetExpansions.ck_output_kernel_complete(output_key)

  for first_sideband in -4:6, second_sideband in -4:6
    sidebands = [first_sideband, second_sideband]
    endpoint_value = FloquetExpansions.ck_endpoint_ordered_sideband_coefficient(
      endpoint, [1, 2], sidebands; inverse_weight=ck_output_inverse_weight
    )
    output_value = FloquetExpansions.ck_output_ordered_sideband_coefficient(
      output, [1, 2], sidebands; inverse_weight=ck_output_inverse_weight
    )
    @test output_value == endpoint_value
  end
end

@testset "drift-only homological cuts become scalar output constraints" begin
  zero_component = zero(CKOutputExact)
  value = one(CKOutputExact)
  vertices = [FloquetExpansions.ck_drift_vertex(2), FloquetExpansions.ck_jump_vertex(1, -1)]
  endpoint_key = FloquetExpansions.CKEndpointKey(
    vertices,
    [FloquetExpansions.CKEndpointConstraint(0, 2)],
    [FloquetExpansions.CKEndpointConstraint(0, 1)],
  )
  endpoint = FloquetExpansions.CKEndpointKernel(Dict(endpoint_key => value), zero_component)
  output = FloquetExpansions.ck_output_kernel(endpoint)
  output_key = only(keys(output.terms))

  @test output_key.blocks == [FloquetExpansions.CKOutputBlock(0, 2, 0, 1)]
  @test output_key.model_constraints == [FloquetExpansions.CKOutputConstraint(1, 1, 1)]
  @test output_key.resolvent_constraints == [FloquetExpansions.CKOutputConstraint(1, 0, 2)]
  @test FloquetExpansions.ck_output_kernel_complete(output_key)

  for sideband in -4:4
    endpoint_value = FloquetExpansions.ck_endpoint_ordered_sideband_coefficient(
      endpoint, [1], [sideband]; inverse_weight=ck_output_inverse_weight
    )
    output_value = FloquetExpansions.ck_output_ordered_sideband_coefficient(
      output, [1], [sideband]; inverse_weight=ck_output_inverse_weight
    )
    @test output_value == endpoint_value
  end
end

@testset "endpoint composition blocks remain distinct in output space" begin
  zero_component = zero(CKOutputExact)
  endpoint_key = FloquetExpansions.CKEndpointKey(
    [FloquetExpansions.ck_jump_vertex(1, 1), FloquetExpansions.ck_jump_vertex(2, 2)],
    [
      FloquetExpansions.CKEndpointConstraint(0, 1),
      FloquetExpansions.CKEndpointConstraint(1, 2),
    ],
    FloquetExpansions.CKEndpointConstraint[],
  )
  endpoint = FloquetExpansions.CKEndpointKernel(
    Dict(endpoint_key => one(CKOutputExact)), zero_component
  )
  output = FloquetExpansions.ck_output_kernel(endpoint)
  output_key = only(keys(output.terms))

  @test output_key.blocks == [
    FloquetExpansions.CKOutputBlock(0, 1, 0, 1), FloquetExpansions.CKOutputBlock(1, 2, 1, 2)
  ]
  @test output_key.model_constraints == [
    FloquetExpansions.CKOutputConstraint(1, 1, 1),
    FloquetExpansions.CKOutputConstraint(2, 2, 2),
  ]
  @test FloquetExpansions.ck_output_kernel_complete(output_key)

  @test FloquetExpansions.ck_output_ordered_sideband_coefficient(
    output, [1, 2], [1, 2]; inverse_weight=ck_output_inverse_weight
  ) == one(CKOutputExact)
  @test iszero(
    FloquetExpansions.ck_output_ordered_sideband_coefficient(
      output, [1, 2], [0, 3]; inverse_weight=ck_output_inverse_weight
    ),
  )
end

@testset "overlapping endpoint block starts are rejected" begin
  endpoint_key = FloquetExpansions.CKEndpointKey(
    [
      FloquetExpansions.ck_jump_vertex(1, 0),
      FloquetExpansions.ck_jump_vertex(1, 0),
      FloquetExpansions.ck_jump_vertex(1, 0),
    ],
    [
      FloquetExpansions.CKEndpointConstraint(0, 2),
      FloquetExpansions.CKEndpointConstraint(1, 3),
    ],
    FloquetExpansions.CKEndpointConstraint[],
  )
  @test_throws ArgumentError FloquetExpansions.ck_output_kernel_key(endpoint_key)
end

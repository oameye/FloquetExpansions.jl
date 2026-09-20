using Test
using FloquetExpansions

const CKOutputRealizationExact = Complex{Rational{Int}}
const ck_output_realization_im = CKOutputRealizationExact(0 // 1, 1 // 1)

function ck_output_realization_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_output_realization_im / mismatch
end

function ck_output_realization_direct_weight(key, sidebands)
  state = FloquetExpansions.CKOutputKernel(
    Dict(key => one(CKOutputRealizationExact)), zero(CKOutputRealizationExact)
  )
  return FloquetExpansions.ck_output_ordered_sideband_coefficient(
    state,
    key.output_channels,
    sidebands;
    inverse_weight=ck_output_realization_inverse_weight,
  )
end

function ck_output_realization_time_weight(key, sidebands)
  realization = FloquetExpansions.ck_output_time_realization(key, ck_output_realization_im)
  return FloquetExpansions.ck_output_time_fourier_weight(
    realization, sidebands, ck_output_realization_im
  )
end

@testset "finite output realization reproduces free cumulative-sideband constraints" begin
  key = FloquetExpansions.CKOutputKernelKey(
    [1, 2],
    [FloquetExpansions.CKOutputBlock(0, 2, 0, 2)],
    FloquetExpansions.CKOutputConstraint{Int}[],
    [
      FloquetExpansions.CKOutputConstraint(1, 1, 2),
      FloquetExpansions.CKOutputConstraint(1, 1, 5),
      FloquetExpansions.CKOutputConstraint(1, 2, -1),
    ],
  )
  realization = FloquetExpansions.ck_output_time_realization(key, ck_output_realization_im)

  @test realization.valid
  @test realization.scalar == one(CKOutputRealizationExact)
  @test length(realization.coordinates) == 2

  for first_sideband in -4:7, second_sideband in -4:7
    sidebands = [first_sideband, second_sideband]
    @test ck_output_realization_time_weight(key, sidebands) ==
      ck_output_realization_direct_weight(key, sidebands)
  end
end

@testset "model-fixed cumulative coordinate becomes a finite physical phase" begin
  key = FloquetExpansions.CKOutputKernelKey(
    [1],
    [FloquetExpansions.CKOutputBlock(0, 1, 0, 1)],
    [FloquetExpansions.CKOutputConstraint(1, 1, 2)],
    [FloquetExpansions.CKOutputConstraint(1, 1, 5)],
  )
  realization = FloquetExpansions.ck_output_time_realization(key, ck_output_realization_im)

  @test realization.valid
  @test realization.scalar == ck_output_realization_im / 3
  @test realization.coordinates[1].kernel.phases ==
    [FloquetExpansions.CKPhaseComponent(2, one(CKOutputRealizationExact))]
  @test isempty(realization.coordinates[1].kernel.homological)

  for sideband in -4:8
    sidebands = [sideband]
    @test ck_output_realization_time_weight(key, sidebands) ==
      ck_output_realization_direct_weight(key, sidebands)
  end
end

@testset "model-resolvent collisions and inconsistent model cuts are invalid" begin
  collision = FloquetExpansions.CKOutputKernelKey(
    [1],
    [FloquetExpansions.CKOutputBlock(0, 1, 0, 1)],
    [FloquetExpansions.CKOutputConstraint(1, 1, 5)],
    [FloquetExpansions.CKOutputConstraint(1, 1, 5)],
  )
  collision_realization =
    FloquetExpansions.ck_output_time_realization(collision, ck_output_realization_im)
  @test !collision_realization.valid

  inconsistent = FloquetExpansions.CKOutputKernelKey(
    [1],
    [FloquetExpansions.CKOutputBlock(0, 1, 0, 1)],
    [
      FloquetExpansions.CKOutputConstraint(1, 1, 2),
      FloquetExpansions.CKOutputConstraint(1, 1, 3),
    ],
    FloquetExpansions.CKOutputConstraint{Int}[],
  )
  inconsistent_realization =
    FloquetExpansions.ck_output_time_realization(inconsistent, ck_output_realization_im)
  @test !inconsistent_realization.valid

  for sideband in -4:8
    @test ck_output_realization_time_weight(collision, [sideband]) ==
      ck_output_realization_direct_weight(collision, [sideband]) ==
      zero(CKOutputRealizationExact)
    @test ck_output_realization_time_weight(inconsistent, [sideband]) ==
      ck_output_realization_direct_weight(inconsistent, [sideband]) ==
      zero(CKOutputRealizationExact)
  end
end

@testset "drift-only scalar constraints are reduced before physical output coordinates" begin
  scalar_resolvent = FloquetExpansions.CKOutputKernelKey(
    [1],
    [FloquetExpansions.CKOutputBlock(0, 1, 0, 1)],
    [FloquetExpansions.CKOutputConstraint(1, 1, 0)],
    [FloquetExpansions.CKOutputConstraint(1, 0, 2)],
  )
  realization =
    FloquetExpansions.ck_output_time_realization(scalar_resolvent, ck_output_realization_im)
  @test realization.valid
  @test realization.scalar == ck_output_realization_im / 2

  invalid_scalar_model = FloquetExpansions.CKOutputKernelKey(
    [1],
    [FloquetExpansions.CKOutputBlock(0, 1, 0, 1)],
    [
      FloquetExpansions.CKOutputConstraint(1, 0, 1),
      FloquetExpansions.CKOutputConstraint(1, 1, 0),
    ],
    FloquetExpansions.CKOutputConstraint{Int}[],
  )
  invalid_realization = FloquetExpansions.ck_output_time_realization(
    invalid_scalar_model, ck_output_realization_im
  )
  @test !invalid_realization.valid

  for sideband in -4:4
    @test ck_output_realization_time_weight(scalar_resolvent, [sideband]) ==
      ck_output_realization_direct_weight(scalar_resolvent, [sideband])
    @test ck_output_realization_time_weight(invalid_scalar_model, [sideband]) ==
      ck_output_realization_direct_weight(invalid_scalar_model, [sideband]) ==
      zero(CKOutputRealizationExact)
  end
end

@testset "cumulative output sidebands reset at endpoint-composition block boundaries" begin
  key = FloquetExpansions.CKOutputKernelKey(
    [1, 2],
    [
      FloquetExpansions.CKOutputBlock(0, 1, 0, 1),
      FloquetExpansions.CKOutputBlock(1, 2, 1, 2),
    ],
    FloquetExpansions.CKOutputConstraint{Int}[],
    [
      FloquetExpansions.CKOutputConstraint(1, 1, 2),
      FloquetExpansions.CKOutputConstraint(2, 2, 5),
    ],
  )
  realization = FloquetExpansions.ck_output_time_realization(key, ck_output_realization_im)
  @test realization.valid
  @test FloquetExpansions.ck_output_time_cumulative_sideband(
    realization, realization.coordinates[1], [3, 7]
  ) == 3
  @test FloquetExpansions.ck_output_time_cumulative_sideband(
    realization, realization.coordinates[2], [3, 7]
  ) == 7

  for first_sideband in -4:7, second_sideband in -4:7
    sidebands = [first_sideband, second_sideband]
    @test ck_output_realization_time_weight(key, sidebands) ==
      ck_output_realization_direct_weight(key, sidebands)
  end
end

@testset "incomplete physical output directions are rejected" begin
  key = FloquetExpansions.CKOutputKernelKey(
    [1],
    [FloquetExpansions.CKOutputBlock(0, 1, 0, 1)],
    FloquetExpansions.CKOutputConstraint{Int}[],
    FloquetExpansions.CKOutputConstraint{Int}[],
  )
  @test !FloquetExpansions.ck_output_kernel_complete(key)
  @test !FloquetExpansions.ck_output_time_realization(key, ck_output_realization_im).valid
end

using Test

const CKKernelBridgeExact = Complex{Rational{Int}}
const ck_kernel_bridge_im = CKKernelBridgeExact(0 // 1, 1 // 1)

struct PhysicalTwoJumpKernel
  first_harmonic::Int
  second_harmonic::Int
end

function two_jump_period_polynomial(first_frequency::Int, second_frequency::Int)
  if iszero(first_frequency) && iszero(second_frequency)
    return CKKernelBridgeExact[0, 0, 1 // 2]
  end

  if iszero(first_frequency)
    return CKKernelBridgeExact[0, ck_kernel_bridge_im / second_frequency]
  end

  first_integral =
    iszero(second_frequency) ? one(CKKernelBridgeExact) : zero(CKKernelBridgeExact)
  second_integral = if iszero(first_frequency + second_frequency)
    one(CKKernelBridgeExact)
  else
    zero(CKKernelBridgeExact)
  end
  linear = (first_integral - second_integral) / (ck_kernel_bridge_im * first_frequency)
  return CKKernelBridgeExact[0, linear]
end

function two_jump_sideband_period(
  first_harmonic::Int, second_harmonic::Int, first_sideband::Int, second_sideband::Int
)
  first_frequency = first_harmonic - first_sideband
  second_frequency = second_harmonic - second_sideband
  return two_jump_period_polynomial(first_frequency, second_frequency)
end

function two_jump_kernel_gram(left::PhysicalTwoJumpKernel, right::PhysicalTwoJumpKernel)
  return two_jump_period_polynomial(
    left.first_harmonic - right.first_harmonic, left.second_harmonic - right.second_harmonic
  )
end

@testset "two-jump physical kernel has unbounded intermediate output-sideband support" begin
  first_harmonic = 1
  second_harmonic = -1
  total_harmonic = first_harmonic + second_harmonic

  for first_sideband in -32:32
    second_sideband = total_harmonic - first_sideband
    polynomial = two_jump_sideband_period(
      first_harmonic, second_harmonic, first_sideband, second_sideband
    )
    @test !all(iszero, polynomial)

    intermediate_mismatch = first_harmonic - first_sideband
    if iszero(intermediate_mismatch)
      @test polynomial == CKKernelBridgeExact[0, 0, 1 // 2]
    else
      @test polynomial ==
        CKKernelBridgeExact[0, ck_kernel_bridge_im / intermediate_mismatch]
    end
  end

  for sideband_cutoff in (1, 2, 4, 8, 16, 32)
    first_sideband = first_harmonic + sideband_cutoff + 1
    second_sideband = total_harmonic - first_sideband
    polynomial = two_jump_sideband_period(
      first_harmonic, second_harmonic, first_sideband, second_sideband
    )
    @test abs(first_sideband) > sideband_cutoff
    @test !all(iszero, polynomial)
    @test polynomial[2] == -ck_kernel_bridge_im / (sideband_cutoff + 1)
  end
end

@testset "finite time-domain simplex kernels quotient the infinite sideband expansion" begin
  support = -1:1
  kernels = [
    PhysicalTwoJumpKernel(first_harmonic, second_harmonic) for first_harmonic in support for
    second_harmonic in support
  ]
  @test length(kernels) == 9

  @test two_jump_kernel_gram(PhysicalTwoJumpKernel(1, -1), PhysicalTwoJumpKernel(0, 0)) ==
    CKKernelBridgeExact[0, ck_kernel_bridge_im]
  @test two_jump_kernel_gram(PhysicalTwoJumpKernel(1, -1), PhysicalTwoJumpKernel(-1, 1)) ==
    CKKernelBridgeExact[0, ck_kernel_bridge_im / 2]
  @test two_jump_kernel_gram(PhysicalTwoJumpKernel(-1, 1), PhysicalTwoJumpKernel(1, -1)) ==
    CKKernelBridgeExact[0, -ck_kernel_bridge_im / 2]

  diagonal = two_jump_kernel_gram(
    PhysicalTwoJumpKernel(1, -1), PhysicalTwoJumpKernel(1, -1)
  )
  @test diagonal == CKKernelBridgeExact[0, 0, 1 // 2]
end

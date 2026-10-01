using Test
using LinearAlgebra: I, eigvals, Hermitian, norm
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions

@testset "native static slot births channels from a PSD dark residual" begin
  residual = ComplexF64[2 1im 0; -1im 1 0; 0 0 0]
  known = zeros(ComplexF64, 3, 3)
  solution = @inferred FE.native_static_solve(
    residual, known, zeros(ComplexF64, 3, 0), Matrix{ComplexF64}[]
  )
  @test solution.canonical
  @test size(solution.correction) == (3, 0)
  @test size(solution.newborn, 2) == 2
  @test solution.newborn * solution.newborn' ≈ residual atol = 1e-10
  @test solution.coefficient ≈ residual atol = 1e-10
end

@testset "native static slot lifts a bright residual onto active amplitudes" begin
  active = reshape(ComplexF64[1, 0, 0], :, 1)
  correction = reshape(ComplexF64[0.3, 0.2 - 0.1im, -0.4im], :, 1)
  known = ComplexF64[0 0 0; 0 0.5 0; 0 0 0]
  residual = known + active * correction' + correction * active'
  solution = FE.native_static_solve(residual, known, active, Matrix{ComplexF64}[])
  @test size(solution.newborn, 2) == 0
  @test norm(solution.dark_residual) <= 1e-12
  reconstructed = known + active * solution.correction' + solution.correction * active'
  @test reconstructed ≈ residual atol = 1e-10
  @test solution.coefficient ≈ residual atol = 1e-10
end

@testset "native static slot uses the full affine PSD slice" begin
  # The minimum-norm point of delta + t*I is indefinite; t = 2 reaches the PSD cone.
  residual = ComplexF64[-2 0; 0 1]
  known = zeros(ComplexF64, 2, 2)
  gauge = [Matrix{ComplexF64}(I, 2, 2)]
  solution = FE.native_static_solve(residual, known, zeros(ComplexF64, 2, 0), gauge)
  @test !solution.canonical
  @test solution.iterations > 0
  @test solution.coordinates ≈ [2.0] atol = 1e-8
  @test minimum(eigvals(Hermitian(solution.dark_residual))) >= -1e-8
  @test solution.coefficient ≈ residual + 2 * gauge[1] atol = 1e-7
  @test solution.newborn * solution.newborn' ≈ ComplexF64[0 0; 0 3] atol = 1e-7

  @test_throws ArgumentError FE.native_static_solve(
    residual, known, zeros(ComplexF64, 2, 0), Matrix{ComplexF64}[]
  )
  @test_throws DimensionMismatch FE.native_static_solve(
    residual, zeros(ComplexF64, 3, 3), zeros(ComplexF64, 2, 0), gauge
  )
end

using Test
using FloquetExpansions
using LinearAlgebra

@testset "shared Gram tangent kernel" begin
  B = ComplexF64[2 0; 0 1; 0 0]
  Dseed = ComplexF64[0.3 + 0.2im -0.1im; 0.4 0.2 - 0.3im; 0.5im -0.7]
  target = B * Dseed' + Dseed * B'

  D = FloquetExpansions.gram_tangent_lift(B, target)
  @test B * D' + D * B' ≈ target atol = 1e-11 rtol = 1e-11

  frame = FloquetExpansions.gram_active_frame(B)
  transformed = frame.active' * D * frame.channel
  Σ = Diagonal(frame.singular)
  @test Σ * transformed ≈ transformed' * Σ atol = 1e-11 rtol = 1e-11

  # Redundant channel columns are reduced to the active right-singular subspace without changing
  # the tangent Gram coefficient.
  redundant = hcat(B, B[:, 1] + 2B[:, 2])
  Dseed_redundant = hcat(Dseed, Dseed[:, 1] - Dseed[:, 2])
  redundant_target = redundant * Dseed_redundant' + Dseed_redundant * redundant'
  Dreduced = FloquetExpansions.gram_tangent_lift(redundant, redundant_target)
  @test redundant * Dreduced' + Dreduced * redundant' ≈ redundant_target atol = 1e-10 rtol =
    1e-10
  @test FloquetExpansions.gram_active_frame(redundant).rank == 2

  # A dark-dark coefficient cannot be generated linearly from the existing active columns.
  obstruction = copy(target)
  obstruction[3, 3] += 1
  @test_throws ArgumentError FloquetExpansions.gram_tangent_lift(B, obstruction)
end

@testset "shared PSD Gram factor" begin
  A = ComplexF64[1 + 0.2im 0.3; -0.4im 0.7 - 0.1im; 0.2 0.5im]
  form = A * A'
  factor = FloquetExpansions.positive_gram_factor(form)
  @test factor * factor' ≈ form atol = 1e-11 rtol = 1e-11

  zero_factor = FloquetExpansions.positive_gram_factor(zeros(ComplexF64, 3, 3))
  @test size(zero_factor) == (3, 0)

  indefinite = ComplexF64[1 0; 0 -0.1]
  @test_throws ArgumentError FloquetExpansions.positive_gram_factor(indefinite)
end

@testset "shared affine PSD section" begin
  # Preserve the canonical minimum-norm affine section whenever it is already PSD.
  canonical_delta = ComplexF64[1 0; 0 2]
  canonical_phi = reshape(Float64[1, 0, 0, 0], 4, 1)
  canonical = FloquetExpansions.positive_affine_section(canonical_delta, canonical_phi)
  @test canonical.canonical
  @test canonical.iterations == 0
  @test canonical.coordinates ≈ [-1.0] atol = 1e-12 rtol = 1e-12
  @test canonical.form ≈ ComplexF64[0 0; 0 2] atol = 1e-12 rtol = 1e-12

  # The minimum-norm affine representative is indefinite, but the same affine line reaches the
  # PSD cone at diag(0, 3). This is the branch the old pseudoinverse-only native solve missed.
  feasible_delta = ComplexF64[-2 0; 0 1]
  feasible_phi = reshape(Float64[1, 1, 0, 0], 4, 1)
  feasible = FloquetExpansions.positive_affine_section(feasible_delta, feasible_phi)
  @test !feasible.canonical
  @test feasible.iterations > 0
  @test feasible.coordinates ≈ [2.0] atol = 1e-8 rtol = 1e-8
  @test feasible.form ≈ ComplexF64[0 0; 0 3] atol = 1e-8 rtol = 1e-8
  @test minimum(eigvals(Hermitian(feasible.form))) >= -1e-9

  # No displacement of the second diagonal entry can repair a fixed negative first eigenvalue.
  infeasible_delta = ComplexF64[-1 0; 0 1]
  infeasible_phi = reshape(Float64[0, 1, 0, 0], 4, 1)
  @test_throws ArgumentError FloquetExpansions.positive_affine_section(
    infeasible_delta, infeasible_phi; maxiter=200
  )
end

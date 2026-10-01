using Test
using LinearAlgebra: I, norm, tr
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions

@testset "dense Liouville representation round-trips GKSL coordinates" begin
  d = 3
  rep = FE.DenseLiouvilleRepresentation(d)
  @test length(rep.basis) == d^2 - 1
  @test all(abs(tr(F)) <= 1e-14 for F in rep.basis)
  gram = [tr(F' * G) for F in rep.basis, G in rep.basis]
  @test gram ≈ Matrix{ComplexF64}(I, d^2 - 1, d^2 - 1) atol = 1e-12

  H = ComplexF64[0.3 0.1im 0; -0.1im -0.2 0.4; 0 0.4 0.5]
  H -= tr(H) / d * I
  X = ComplexF64[(j + 2k) / 7 + im * (j - k) / 5 for j in 1:(d ^ 2 - 1), k in 1:2]
  C = X * X'
  L = FE.native_gksl(rep, H, C)
  @test FE.native_kossakowski(rep, L) ≈ C atol = 1e-12
  @test FE.native_hamiltonian(rep, L) ≈ H atol = 1e-12
  @test norm(L' * vec(Matrix{ComplexF64}(I, d, d))) <= 1e-12
end

@testset "dense gauge directions are Hermiticity and trace preserving" begin
  d = 2
  rep = FE.DenseLiouvilleRepresentation(d)
  directions = FE.native_gauge_directions(rep)
  @test length(directions) == d^4 - d^2
  identity = vec(Matrix{ComplexF64}(I, d, d))
  @test all(norm(G' * identity) <= 1e-12 for G in directions)
  ρ = ComplexF64[0.6 0.2im; -0.2im 0.4]
  @test all(
    norm(reshape(G * vec(ρ), d, d) - reshape(G * vec(ρ), d, d)') <= 1e-12 for
    G in directions
  )
  @test_throws ArgumentError FE.DenseLiouvilleRepresentation(1)
end

@testset "generic native static step on the dense representation" begin
  d = 2
  rep = FE.DenseLiouvilleRepresentation(d)
  σm = ComplexF64[0 0; 1 0]
  L0 = FE.native_gksl(rep, ComplexF64[0.5 0; 0 -0.5], zeros(ComplexF64, 3, 3))
  L0 +=
    FE.dense_sandwich(σm, σm) - FE.dense_left(σm' * σm) / 2 - FE.dense_right(σm' * σm) / 2
  active = reshape(ComplexF64[tr(F' * σm) for F in rep.basis], :, 1)

  # A pure bright perturbation needs no gauge and no newborn channel.
  correction = 0.1 * ComplexF64[1, -1im, 0.5]
  Cbright = active * correction' + correction * active'
  residual = FE.native_gksl(rep, zeros(ComplexF64, d, d), Cbright)
  step = FE.native_static_step(
    rep, FE.NoHomologicalInverse(), L0, residual, zeros(ComplexF64, 3, 3), active, 1e-8
  )
  @test step isa FE.NativeStaticStep{Matrix{ComplexF64},Matrix{ComplexF64}}
  @test size(step.solution.newborn, 2) == 0
  @test FE.native_kossakowski(rep, step.E) ≈ step.solution.coefficient atol = 1e-9
  @test step.E ≈ FE.native_gksl(rep, step.H, step.solution.coefficient) atol = 1e-8
end

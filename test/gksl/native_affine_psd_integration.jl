using Test
using LinearAlgebra
using FloquetExpansions

@testset "affine PSD section preserves canonical branch" begin
  delta = ComplexF64[1 0; 0 2]
  phi = reshape(Float64[1, 0, 0, 0], 4, 1)
  section = FloquetExpansions.positive_affine_section(delta, phi)

  @test section.canonical
  @test section.coordinates == [-1.0]
  @test section.form == ComplexF64[0 0; 0 2]
end

@testset "affine PSD section rejects disjoint slice" begin
  delta = ComplexF64[-1 0; 0 1]
  phi = reshape(Float64[0, 1, 0, 0], 4, 1)
  @test_throws ArgumentError FloquetExpansions.positive_affine_section(
    delta, phi; maxiter=200
  )
end

using Test
using LinearAlgebra
using FloquetExpansions

# Focused research fixture for #368. The full BF/HD recurrence still lives in the frozen stacked
# tests; this file certifies the new affine-cone section independently before it is wired into the
# native static step.

@testset "native GKSL affine static section geometry" begin
  delta = ComplexF64[-2 0; 0 1]
  phi = reshape(Float64[1, 1, 0, 0], 4, 1)

  canonical_coordinate = -pinv(phi) * Float64[-2, 1, 0, 0]
  canonical_residual = delta + canonical_coordinate[1] * Matrix{ComplexF64}(I, 2, 2)
  @test minimum(eigvals(Hermitian(canonical_residual))) < 0

  section = FloquetExpansions.positive_affine_section(delta, phi)
  @test !section.canonical
  @test section.coordinates ≈ [2.0] atol = 1e-8 rtol = 1e-8
  @test section.form ≈ ComplexF64[0 0; 0 3] atol = 1e-8 rtol = 1e-8

  newborn = FloquetExpansions.positive_gram_factor(section.form)
  @test newborn * newborn' ≈ section.form atol = 1e-8 rtol = 1e-8
end

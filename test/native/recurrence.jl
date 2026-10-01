using Test
using LinearAlgebra: norm, tr
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions

function driven_qubit_harmonics(rep)
  T = ComplexF64
  σx = T[0 1; 1 0]
  σz = T[1 0; 0 -1]
  σm = T[0 0; 1 0]
  zero3 = zeros(T, 3, 3)
  coefficient(X) = [tr(F' * X) for F in rep.basis]
  loss = 0.6 * coefficient(σm) * coefficient(σm)'
  L = Dict{Int,Matrix{T}}(
    0 => FE.native_gksl(rep, 0.3 * σz, loss),
    1 => FE.native_gksl(rep, 0.4 * σx, zero3),
    -1 => FE.native_gksl(rep, 0.4 * σx, zero3),
  )
  return L, reshape(sqrt(0.6) * coefficient(σm), :, 1)
end

@testset "native recurrence: BF and HD agree through order 4" begin
  rep = FE.DenseLiouvilleRepresentation(2)
  L, leading = driven_qubit_harmonics(rep)
  inverse = FE.NoHomologicalInverse()
  bf = @inferred FE.native_recurrence(FE.BlochFeshbach(), rep, inverse, L, leading, 4, 1e-8)
  hd = FE.native_recurrence(FE.HoriDeprit(), rep, inverse, L, leading, 4, 1e-8)

  @test bf isa FE.NativeRecurrence{Matrix{ComplexF64},Matrix{ComplexF64}}
  @test length(bf.E) == 5
  @test all(norm(bf.E[k] - hd.E[k]) <= 1e-7 * max(1.0, norm(bf.E[k])) for k in 1:5)
  @test all(norm(bf.H[k] - hd.H[k]) <= 1e-7 * max(1.0, norm(bf.H[k])) for k in 1:5)
  for order in 0:4
    @test FE.native_kossakowski(rep, bf.E[order + 1]) ≈
      FE.gram_coefficient(bf.channels, order, 3) atol = 1e-7
  end

  # Truncating the retained order does not change the lower coefficients.
  bf2 = FE.native_recurrence(FE.BlochFeshbach(), rep, inverse, L, leading, 2, 1e-8)
  @test all(norm(bf2.E[k] - bf.E[k]) <= 1e-9 * max(1.0, norm(bf.E[k])) for k in 1:3)
end

@testset "native recurrence: order zero and input validation" begin
  rep = FE.DenseLiouvilleRepresentation(2)
  L, leading = driven_qubit_harmonics(rep)
  inverse = FE.NoHomologicalInverse()
  zeroth = FE.native_recurrence(FE.BlochFeshbach(), rep, inverse, L, leading, 0, 1e-8)
  @test zeroth.E == [L[0]]
  @test isempty(zeroth.S)

  @test_throws ArgumentError FE.native_recurrence(
    FE.BlochFeshbach(), rep, inverse, L, leading, -1, 1e-8
  )
  without_average = Dict(k => v for (k, v) in L if k != 0)
  @test_throws ArgumentError FE.native_recurrence(
    FE.HoriDeprit(), rep, inverse, without_average, leading, 2, 1e-8
  )
  # Leading columns that do not reproduce the averaged Kossakowski form are rejected.
  @test_throws ArgumentError FE.native_recurrence(
    FE.BlochFeshbach(), rep, inverse, L, 2 * leading, 2, 1e-8
  )
end

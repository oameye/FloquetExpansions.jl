include("native_bf_prototype.jl")

@testset "native GKLS static step uses full affine PSD slice" begin
  d = 3
  H0 = ComplexF64[0 0 0; 0 0.7 0; 0 0 1.9]
  R = zeros(ComplexF64, d, d)
  R[1, 2] = 0.2
  L0 = native_lindblad(H0, [R])
  active = reshape(native_tlcoef(R, d), :, 1)
  _, dark_basis = native_active_split(active)

  # This quotient defect is deliberately chosen so that the Moore--Penrose section is indefinite,
  # although the same physical static-gauge affine slice intersects the PSD cone.
  quotient_delta = Diagonal(ComplexF64[1, 1, 0, 0, 0, 0, 0]) |> Matrix
  cdefect = dark_basis * quotient_delta * dark_basis'
  Vhat = native_from_Hc(zeros(ComplexF64, d, d), cdefect, d)
  known = zeros(ComplexF64, d^2 - 1, d^2 - 1)

  gauge_basis = native_gauge_algebra(d)
  Φ = native_phi_matrix(L0, d, dark_basis, gauge_basis)
  canonical_coordinates = -pinv(Φ; rtol=1e-10) * native_hvec(quotient_delta)
  canonical_residual = native_hmat(
    native_hvec(quotient_delta) + Φ * canonical_coordinates, size(quotient_delta, 1)
  )
  @test minimum(eigvals(Hermitian(canonical_residual))) < -0.3

  dark = native_dark_solve(L0, Vhat, known, active, d)
  @test !dark.canonical
  @test dark.iterations > 0
  @test minimum(eigvals(Hermitian(dark.dark_residual))) >= -1e-8

  step = native_static_step(L0, Vhat, known, active, d)
  reconstructed =
    step.correction * active' + active * step.correction' + step.newborn * step.newborn'
  @test native_kossakowski(step.E, d) ≈ reconstructed atol = 1e-8 rtol = 1e-8
end

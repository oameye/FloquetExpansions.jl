include("native_hd_arbitrary.jl")

@testset "native GKSL affine PSD branch: periodic qutrit BF equals HD" begin
  d = 3
  H0 = Diagonal(ComplexF64[0, 0.7, 1.9]) |> Matrix
  H1 = ComplexF64[
    -0.0069425955001163-0.0874101272620143im 0.0735434729136112+0.0919593343134481im 0.07894681944654-0.191530360055452im
    -0.0413283199704676+0.098827455433883im -0.0596846442595027-0.0750679764310141im -0.0137727340088244-0.0655127946832933im
    -0.0726542411220522-0.1621016570238948im -0.0713207301005824-0.0173090542618764im -0.0340050450724749-0.0297193811120841im
  ]
  R0 = zeros(ComplexF64, d, d)
  R0[1, 2] = 0.2
  R1 = zeros(ComplexF64, d, d)
  R1[2, 3] = 0.08 + 0.024im
  model = NativeModel(
    d, NativeFS(0 => H0, 1 => H1, -1 => H1'), [NativeFS(0 => R0, 1 => R1)]
  )

  bf = native_bf_arbitrary(model, 2)
  hd = native_hd_arbitrary(model, 2)
  native_compare_arbitrary_states(bf, hd, 2; tol=2e-6)
  @test maximum(native_hd_equation_defects(hd, model, 2)) <= 2e-6
  @test maximum(native_hd_intrinsic_defects(hd, 2, d^2)) <= 2e-6

  # Reconstruct the order-two retained static problem before the accepted static slot. The
  # Moore--Penrose section is genuinely indefinite, while the full affine slice is feasible.
  L = native_liouvillian_harmonics(model)
  L0 = native_fsavg(L, d^2)
  offset = native_bf_intrinsic_offset(bf.Y, bf.S[1:1], 2, d^2)
  Ynbase = native_set_static(bf.Y[3], -offset)
  Vhat = native_fsavg(native_bf_static_residual(L, bf.Y, bf.E, 2, Ynbase), d^2)
  known = native_known_gram(bf.channels, 2, d^2 - 1)
  _, active = native_active_channels(bf.channels, 2, d^2 - 1)
  _, dark_basis = native_active_split(active)
  delta = native_hermitian(dark_basis' * (native_kossakowski(Vhat, d) - known) * dark_basis)
  Φ = native_phi_matrix(L0, d, dark_basis, native_gauge_algebra(d))
  canonical_coordinates = -pinv(Φ; rtol=1e-10) * native_hvec(delta)
  canonical_residual = native_hmat(
    native_hvec(delta) + Φ * canonical_coordinates, size(delta, 1)
  )
  @test minimum(eigvals(Hermitian(canonical_residual))) < -1e-6

  resolved = native_dark_solve(L0, Vhat, known, active, d)
  @test !resolved.canonical
  @test resolved.iterations > 0
  @test minimum(eigvals(Hermitian(resolved.dark_residual))) >= -1e-8

  ε = 1e-3
  finite_bf = native_bf_finite(bf, d, ε)
  finite_hd = native_hd_finite(hd, d, ε)
  @test norm(finite_bf - finite_hd) <= 2e-6
  @test minimum(eigvals(Hermitian(native_hermitian(native_kossakowski(finite_bf, d))))) >=
    -1e-10
end

include("native_hd_arbitrary.jl")

function native_fock_annihilation(d)
  a = zeros(ComplexF64, d, d)
  for n in 1:(d - 1)
    a[n, n + 1] = sqrt(n)
  end
  return a
end

@testset "native GKLS BF/HD finite-Fock driven Kerr" begin
  d = 3
  Δ = 1 / 2
  χ = 3 / 10
  κ = 4 / 5
  F = 2 / 5 + im / 5

  a = native_fock_annihilation(d)
  n = a' * a
  x = a + a'
  H0 = Δ * n + χ * n * n
  H1 = F * x
  model = NativeModel(
    d, NativeFS(0 => H0, 1 => H1, -1 => H1'), [NativeFS(0 => sqrt(κ) * a)]
  )

  bf2 = native_bf_arbitrary(model, 2)
  bf4 = native_bf_arbitrary(model, 4)
  hd4 = native_hd_arbitrary(model, 4)

  native_compare_arbitrary_states(bf4, hd4, 4; tol=2e-6)
  @test maximum(native_hd_equation_defects(hd4, model, 4)) <= 2e-6
  @test maximum(native_hd_intrinsic_defects(hd4, 4, d^2)) <= 2e-6

  @test all(isapprox(bf2.E[k], bf4.E[k]; atol=2e-7, rtol=2e-7) for k in 1:3)
  @test all(isapprox(bf2.S[k], bf4.S[k]; atol=2e-7, rtol=2e-7) for k in 1:2)
  @test native_channel_prefix_equal(bf2, bf4, 2; atol=2e-7)

  dim = d^2 - 1
  for order in 0:4
    gram = native_gram_coefficient(bf4.channels, order, dim)
    @test isapprox(native_kossakowski(bf4.E[order + 1], d), gram; atol=2e-6, rtol=2e-6)
  end

  ε = 5e-4
  finite_bf = native_bf_finite(bf4, d, ε)
  finite_hd = native_hd_finite(hd4, d, ε)
  @test norm(finite_bf - finite_hd) <= 2e-6
  @test norm(finite_bf' * vec(native_id(d))) <= 2e-8
  cfinite = native_hermitian(native_kossakowski(finite_bf, d))
  @test minimum(eigvals(Hermitian(cfinite))) >= -2e-8
end

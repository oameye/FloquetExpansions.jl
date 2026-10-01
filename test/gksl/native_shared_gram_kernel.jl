using FloquetExpansions

include("native_hd_arbitrary.jl")

# Certify that the research recurrence runs its static slot through the package core: the
# order-four BF/HD fixture is reproduced, and a single step equals a direct call of
# native_static_solve on the same Kossakowski data.
const shared_gram_γ = 0.52
const shared_gram_model = NativeModel(
  2,
  NativeFS(
    0 => 0.21 * σz_native,
    1 => 0.24 * σx_native,
    -1 => 0.24 * σx_native,
    2 => 0.07 * σy_native,
    -2 => 0.07 * σy_native,
  ),
  [
    NativeFS(
      -1 => sqrt(shared_gram_γ) * 0.31 * σx_native,
      0 => sqrt(shared_gram_γ) * 0.43 * σy_native,
      1 => sqrt(shared_gram_γ) * 0.57 * σz_native,
    ),
  ],
)
const shared_gram_bf_baseline = native_bf_arbitrary(shared_gram_model, 4)
const shared_gram_hd_baseline = native_hd_arbitrary(shared_gram_model, 4)

function native_shared_state_equal(a, b, order; tol=2e-7)
  for n in 0:order
    isapprox(a.E[n + 1], b.E[n + 1]; atol=tol, rtol=tol) || return false
    isapprox(a.H[n + 1], b.H[n + 1]; atol=tol, rtol=tol) || return false
  end
  for n in 1:order
    isapprox(a.S[n], b.S[n]; atol=tol, rtol=tol) || return false
  end
  return native_channel_prefix_equal(a, b, order; atol=tol)
end

@testset "native GKLS reuses shared Gram kernel" begin
  bf = native_bf_arbitrary(shared_gram_model, 4)
  hd = native_hd_arbitrary(shared_gram_model, 4)

  @test native_shared_state_equal(bf, shared_gram_bf_baseline, 4)
  @test native_shared_state_equal(hd, shared_gram_hd_baseline, 4)
  native_compare_arbitrary_states(bf, hd, 4; tol=2e-7)

  d = shared_gram_model.d
  L0 = native_fsavg(native_liouvillian_harmonics(shared_gram_model), d^2)
  Vhat = bf.E[3] - (L0 * bf.S[2] - bf.S[2] * L0)
  _, active = native_active_channels(bf.channels, 2, d^2 - 1)
  known = bf.known_gram[2]
  step = native_static_step(L0, Vhat, known, active, d)
  images = NativeCM[native_kossakowski(L0 * G - G * L0, d) for G in native_gauge_algebra(d)]
  direct = FloquetExpansions.native_static_solve(
    native_kossakowski(Vhat, d), known, active, images; tol=1e-8
  )
  @test step.C ≈ direct.coefficient atol = 1e-10 rtol = 1e-10
  @test step.correction ≈ direct.correction atol = 1e-10 rtol = 1e-10
end

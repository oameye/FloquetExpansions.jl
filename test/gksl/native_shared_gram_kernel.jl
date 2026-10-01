using FloquetExpansions

include("native_hd_arbitrary.jl")

# Certify one nontrivial order-four result through the generic shared-kernel adapters, then install
# more-specific NativeCM methods below and require the recurrence to remain coefficientwise
# unchanged. This checks that the package kernel is independent of the adapter dispatch used by
# the research oracle.
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

# These methods are deliberately more specific than the generic adapters in
# native_bf_prototype.jl. Both paths call the same package-level Gram geometry; the specialized
# dispatch makes that production-kernel reuse explicit in the order-four BF/HD fixture.
function native_active_split(A::NativeCM; tol=1e-10)
  frame = FloquetExpansions.gram_active_frame(A; rtol=tol)
  return frame.active, frame.dark
end

function native_tangent_lift(active::NativeCM, target::NativeCM; tol=1e-9)
  return FloquetExpansions.gram_tangent_lift(active, target; rtol=tol)
end

function native_dark_solve(
  L0::NativeCM,
  Vhat::NativeCM,
  known::NativeCM,
  active::NativeCM,
  d::Int;
  tol=1e-9,
  gauge_basis=native_gauge_algebra(d),
)
  frame = FloquetExpansions.gram_active_frame(active; rtol=tol)
  P = frame.dark
  darkdim = size(P, 2)
  if darkdim == 0
    return (
      S=zeros(ComplexF64, d^2, d^2),
      newborn=zeros(ComplexF64, d^2 - 1, 0),
      dark_residual=zeros(ComplexF64, 0, 0),
    )
  end

  delta = native_hermitian(P' * (native_kossakowski(Vhat, d) - known) * P)
  Φ = native_phi_matrix(L0, d, P, gauge_basis)
  coordinates = isempty(gauge_basis) ? Float64[] : -pinv(Φ; rtol=1e-10) * native_hvec(delta)

  S = zeros(ComplexF64, d^2, d^2)
  for i in eachindex(coordinates)
    S += coordinates[i] * gauge_basis[i]
  end

  residual = native_hermitian(delta + P' * native_kossakowski(L0 * S - S * L0, d) * P)
  quotient_factor = FloquetExpansions.positive_gram_factor(residual; rtol=tol)
  newborn = P * quotient_factor
  return (; S, newborn, dark_residual=residual)
end

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

  # The package-kernel adapters are the selected dispatch for the dense native matrices.
  @test String(which(native_tangent_lift, (NativeCM, NativeCM)).file) == @__FILE__
  @test String(which(native_active_split, (NativeCM,)).file) == @__FILE__
end

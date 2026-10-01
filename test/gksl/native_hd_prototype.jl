using Test
using LinearAlgebra

# Reuse only the dense GKLS coordinate/static-solve oracle. The Floquet recurrence below is an
# independent Hori--Deprit construction and is compared against the BF recurrence coefficient by
# coefficient in the same intrinsic total-kick gauge.
include("native_bf_prototype.jl")

function native_fsscale(A::NativeFS, α)
  return NativeFS(k => α * value for (k, value) in A)
end

function native_fscomm(A::NativeFS, B::NativeFS)
  return native_fsadd(native_fsmul(A, B), native_fsmul(B, A), 1.0, -1.0)
end

function native_fsder(A::NativeFS)
  return NativeFS(k => (-im * k) * value for (k, value) in A if k != 0)
end

function native_fs_with_static(A::NativeFS, S::NativeCM)
  result = copy(A)
  result[0] = S
  return result
end

function native_hd_forcing1(L::NativeFS, G1::NativeFS)
  dG1 = native_fsder(G1)
  return native_fsadd(
    native_fsscale(native_fscomm(G1, L), -1.0), native_fsscale(native_fscomm(G1, dG1), 0.5)
  )
end

function native_hd_forcing2(L::NativeFS, G1::NativeFS, G2::NativeFS)
  dG1 = native_fsder(G1)
  dG2 = native_fsder(G2)

  result = native_fsscale(native_fscomm(G2, L), -1.0)
  result = native_fsadd(
    result, native_fsscale(native_fscomm(G1, native_fscomm(G1, L)), 0.5)
  )
  result = native_fsadd(result, native_fsscale(native_fscomm(G1, dG2), 0.5))
  result = native_fsadd(result, native_fsscale(native_fscomm(G2, dG1), 0.5))
  result = native_fsadd(
    result, native_fsscale(native_fscomm(G1, native_fscomm(G1, dG1)), -1.0 / 6.0)
  )
  return result
end

function native_pack_order02(E, S, H, B0, B1, B2, U0, U1, V0, G1osc, G2osc, K2known, d)
  birth0 = NativeCM[]
  for j in axes(B0, 2)
    append!(
      birth0,
      [
        native_operator(B0[:, j], d),
        native_operator(B1[:, j], d),
        native_operator(B2[:, j], d),
      ],
    )
  end

  birth1 = NativeCM[]
  for j in axes(U0, 2)
    append!(birth1, [native_operator(U0[:, j], d), native_operator(U1[:, j], d)])
  end

  birth2 = NativeCM[native_operator(V0[:, j], d) for j in axes(V0, 2)]

  return NativeOrder02(E, S, copy(S), H, birth0, birth1, birth2, G1osc, G2osc, K2known)
end

function native_hd_order02(model::NativeModel; tol=1e-8)
  d = model.d
  nsuper = d^2
  L = native_liouvillian_harmonics(model)
  L0 = native_fsavg(L, nsuper)

  # Order 0: L - ∂τ G1 = E0.
  E0 = L0
  G1osc = native_fsint(L)
  B0 = native_sideband_columns(model)
  c0 = B0 * B0'
  norm(native_kossakowski(E0, d) - c0) <= tol * max(1.0, norm(c0)) || error(
    "HD order-zero sideband columns do not reconstruct the averaged Kossakowski tensor"
  )
  H0 = native_hamiltonian_part(E0, d)

  # Order 1:
  # F1 = -[G1,L] + 1/2 [G1,∂τG1], with the yet-unknown static S1 set to zero
  # when the intrinsic residual is formed. The static solve then supplies S1.
  Vhat1 = native_fsavg(native_hd_forcing1(L, G1osc), nsuper)
  known1 = zeros(ComplexF64, d^2 - 1, d^2 - 1)
  step1 = native_static_step(L0, Vhat1, known1, B0, d; tol)
  S1 = step1.S
  E1 = step1.E
  B1 = step1.correction
  U0 = step1.newborn

  G1 = native_fs_with_static(G1osc, S1)
  forcing1 = native_hd_forcing1(L, G1)
  G2osc = native_fsint(forcing1)

  # Order 2:
  # F2 = -[G2,L] + 1/2[G1,[G1,L]]
  #      + 1/2([G1,∂G2]+[G2,∂G1]) - 1/6[G1,[G1,∂G1]].
  # In the intrinsic total-kick convention <G2> = S2 through this order, so zero static G2
  # directly gives the intrinsic residual Vhat2.
  Vhat2 = native_fsavg(native_hd_forcing2(L, G1, G2osc), nsuper)
  K2known = B1 * B1'
  active2 = hcat(B0, U0)
  step2 = native_static_step(L0, Vhat2, K2known, active2, d; tol)
  S2 = step2.S
  E2 = step2.E

  correction2 = step2.correction
  B2 = correction2[:, 1:size(B0, 2)]
  U1 = correction2[:, (size(B0, 2) + 1):end]
  V0 = step2.newborn

  return native_pack_order02(
    [E0, E1, E2],
    [S1, S2],
    [H0, step1.H, step2.H],
    B0,
    B1,
    B2,
    U0,
    U1,
    V0,
    G1osc,
    G2osc,
    K2known,
    d,
  )
end

function native_compare_order02(bf::NativeOrder02, hd::NativeOrder02; tol=2e-7)
  for n in 1:3
    @test norm(bf.E[n] - hd.E[n]) <= tol * max(1.0, norm(bf.E[n]))
    @test norm(bf.H[n] - hd.H[n]) <= tol * max(1.0, norm(bf.H[n]))
  end
  for n in 1:2
    @test norm(bf.S[n] - hd.S[n]) <= tol * max(1.0, norm(bf.S[n]))
  end
  @test length(bf.birth0) == length(hd.birth0)
  @test length(bf.birth1) == length(hd.birth1)
  @test length(bf.birth2) == length(hd.birth2)
  @test norm(bf.K2known - hd.K2known) <= tol * max(1.0, norm(bf.K2known))
end

@testset "native GKLS HD/BF equality: Hamiltonian reduction through order 2" begin
  model = NativeModel(
    2,
    NativeFS(
      0 => 0.37 * σz_native,
      1 => 0.41 * σx_native,
      -1 => 0.41 * σx_native,
      2 => 0.13 * σy_native,
      -2 => 0.13 * σy_native,
    ),
    NativeFS[],
  )

  bf = native_bf_order02(model)
  hd = native_hd_order02(model)
  native_compare_order02(bf, hd)

  @test norm(hd.S[1]) <= 1e-9
  @test norm(hd.S[2]) <= 1e-9
  @test isempty(hd.birth0)
  @test isempty(hd.birth1)
  @test isempty(hd.birth2)
end

@testset "native GKLS HD/BF equality: driven qubit with static loss" begin
  γ = 0.63
  model = NativeModel(
    2,
    NativeFS(0 => 0.29 * σz_native, 1 => 0.38 * σx_native, -1 => 0.38 * σx_native),
    [NativeFS(0 => sqrt(γ) * σm_native)],
  )

  bf = native_bf_order02(model)
  hd = native_hd_order02(model)
  native_compare_order02(bf, hd)

  ε = 1e-3
  @test norm(native_effective_gkls(bf, ε) - native_effective_gkls(hd, ε)) <= 1e-9
end

@testset "native GKLS HD/BF equality: full-rank physical sidebands" begin
  γ = 0.52
  sideband_jump = NativeFS(
    -1 => sqrt(γ) * 0.31 * σx_native,
    0 => sqrt(γ) * 0.43 * σy_native,
    1 => sqrt(γ) * 0.57 * σz_native,
  )
  model = NativeModel(
    2,
    NativeFS(0 => 0.21 * σz_native, 1 => 0.24 * σx_native, -1 => 0.24 * σx_native),
    [sideband_jump],
  )

  bf = native_bf_order02(model)
  hd = native_hd_order02(model)
  native_compare_order02(bf, hd)

  @test isempty(hd.birth1)
  @test isempty(hd.birth2)
end

@testset "native GKLS static solve: positive order-one channel birth" begin
  d = 2
  L0 = zeros(ComplexF64, d^2, d^2)
  born_jump = 0.37 * σx_native
  Vhat1 = native_lindblad(zeros(ComplexF64, d, d), [born_jump])
  known = zeros(ComplexF64, d^2 - 1, d^2 - 1)
  active = zeros(ComplexF64, d^2 - 1, 0)

  step = native_static_step(L0, Vhat1, known, active, d)

  @test norm(step.S) <= 1e-10
  @test size(step.newborn, 2) == 1
  @test size(step.correction, 2) == 0

  recovered_jump = native_operator(step.newborn[:, 1], d)
  @test norm(step.E - native_lindblad(step.H, [recovered_jump])) <= 1e-9

  # A generator coefficient born at order ε is represented by an amplitude born at order ε^(1/2).
  ε = 1e-4
  finite = native_lindblad(ε * step.H, [sqrt(ε) * recovered_jump])
  @test norm(finite - ε * step.E) <= 1e-10
end

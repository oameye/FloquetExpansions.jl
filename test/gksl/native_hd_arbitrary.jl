include("native_bf_arbitrary.jl")

# Arbitrary-order native Hori--Deprit recurrence, run through the package driver.
#
# The Lie-transform series below is an independent re-implementation of the HD identity. The
# equation-defect checks evaluate the package HD state against it, and the package BF and HD
# drivers share only the static step and the graded channel state.

struct NativeHDState
  E::Vector{NativeCM}
  S::Vector{NativeCM}
  hd_slots::Vector{NativeCM}
  H::Vector{NativeCM}
  G::Vector{NativeFS}
  channels::Vector{NativeChannelSeries}
  known_gram::Vector{NativeCM}
end

function native_hd_fsder(A::NativeFS)
  return NativeFS(k => (-im * k) * value for (k, value) in A if k != 0)
end

function native_hd_series_comm(A, B, N)
  AB = native_arb_series_mul(A, B, N)
  BA = native_arb_series_mul(B, A, N)
  return native_arb_series_add(AB, BA, N; α=1.0, β=-1.0)
end

function native_hd_generator_series(L::NativeFS, G, N)
  Gseries = native_arb_series_zero(N)
  for n in 1:min(N, length(G) - 1)
    Gseries[n + 1] = G[n + 1]
  end

  Lseries = native_arb_series_zero(N)
  Lseries[1] = L

  # exp(-ad_G) L
  effective = [copy(value) for value in Lseries]
  adterm = Lseries
  for k in 1:N
    adterm = native_hd_series_comm(Gseries, adterm, N)
    coefficient = (-1.0)^k / factorial(k)
    effective = native_arb_series_add(
      effective, native_arb_series_scale(adterm, N, coefficient), N
    )
  end

  # Dbar = ε^-1 ∂τ G, so Dbar_n = ∂τ G_{n+1}.
  Dbar = native_arb_series_zero(N)
  for n in 0:N
    index = n + 2
    index <= length(G) || continue
    Dbar[n + 1] = native_hd_fsder(G[index])
  end

  # -(1 - exp(-ad_G))/ad_G Dbar
  dterm = Dbar
  effective = native_arb_series_add(effective, dterm, N; α=1.0, β=-1.0)
  for k in 1:N
    dterm = native_hd_series_comm(Gseries, dterm, N)
    coefficient = -(-1.0)^k / factorial(k + 1)
    effective = native_arb_series_add(
      effective, native_arb_series_scale(dterm, N, coefficient), N
    )
  end
  return effective
end

function native_hd_arbitrary(model::NativeModel, N; tol=1e-8, inverse=NativeNoInverse())
  result = FloquetExpansions.native_recurrence(
    FloquetExpansions.HoriDeprit(),
    native_dense_representation(model.d),
    inverse,
    native_liouvillian_harmonics(model),
    native_sideband_columns(model),
    N,
    tol,
  )
  return NativeHDState(
    result.E,
    result.S,
    result.slots,
    result.H,
    result.kick,
    result.channels,
    result.known_gram,
  )
end

function native_hd_finite(state::NativeHDState, d, ε)
  H = zeros(ComplexF64, d, d)
  for n in eachindex(state.H)
    H += ε^(n - 1) * state.H[n]
  end

  jumps = NativeCM[]
  for channel in state.channels
    R = zeros(ComplexF64, d, d)
    for (k, coefficient) in pairs(channel.coefficients)
      R += ε^(k - 1) * native_operator(coefficient, d)
    end
    push!(jumps, ε^(channel.onset / 2) * R)
  end
  return native_lindblad(H, jumps)
end

function native_hd_static_defect(A::NativeFS, E::NativeCM)
  modes = union(keys(A), (0,))
  defect = 0.0
  for mode in modes
    target = mode == 0 ? E : zero(E)
    defect = max(defect, norm(get(A, mode, zero(E)) - target))
  end
  return defect
end

function native_hd_equation_defects(state::NativeHDState, model::NativeModel, N)
  L = native_liouvillian_harmonics(model)
  defects = Float64[]
  for order in 0:N
    Gcheck = [copy(value) for value in state.G]
    if order == N
      transformed = native_hd_generator_series(L, Gcheck, N)
      push!(Gcheck, native_fsint(transformed[N + 1]))
    end
    transformed = native_hd_generator_series(L, Gcheck, order)
    push!(defects, native_hd_static_defect(transformed[order + 1], state.E[order + 1]))
  end
  return defects
end

function native_hd_intrinsic_defects(state::NativeHDState, N, nsuper)
  Gseries = native_arb_series_zero(N)
  for n in 1:N
    Gseries[n + 1] = state.G[n + 1]
  end
  kick = native_arb_series_exp(Gseries, nsuper, N)
  minusS = native_arb_static_exp(state.S, nsuper, N; sign=-1.0)
  Kosc = native_arb_series_mul(kick, minusS, N)
  logKosc = native_arb_series_log(Kosc, nsuper, N)
  return [norm(native_fsavg(logKosc[n + 1], nsuper)) for n in 1:N]
end

function native_compare_arbitrary_states(bf::NativeBFState, hd::NativeHDState, N; tol=5e-7)
  @test length(bf.E) == N + 1
  @test length(hd.E) == N + 1
  d = isqrt(size(bf.E[1], 1))
  dim = d^2 - 1
  for order in 0:N
    @test norm(bf.E[order + 1] - hd.E[order + 1]) <= tol * max(1.0, norm(bf.E[order + 1]))
    @test norm(bf.H[order + 1] - hd.H[order + 1]) <= tol * max(1.0, norm(bf.H[order + 1]))
    cbf = native_gram_coefficient(bf.channels, order, dim)
    chd = native_gram_coefficient(hd.channels, order, dim)
    @test norm(cbf - chd) <= tol * max(1.0, norm(cbf))
  end
  for order in 1:N
    @test norm(bf.S[order] - hd.S[order]) <= tol * max(1.0, norm(bf.S[order]))
  end
end

@testset "arbitrary native GKLS HD reproduces the frozen low-order oracle" begin
  γ = 0.63
  model = NativeModel(
    2,
    NativeFS(0 => 0.29 * σz_native, 1 => 0.38 * σx_native, -1 => 0.38 * σx_native),
    [NativeFS(0 => sqrt(γ) * σm_native)],
  )
  frozen = native_bf_order02(model)
  hd = native_hd_arbitrary(model, 2)

  @test all(isapprox(hd.E[n], frozen.E[n]; atol=2e-7, rtol=2e-7) for n in 1:3)
  @test all(isapprox(hd.S[n], frozen.S[n]; atol=2e-7, rtol=2e-7) for n in 1:2)
  @test all(isapprox(hd.H[n], frozen.H[n]; atol=2e-7, rtol=2e-7) for n in 1:3)
  @test maximum(native_hd_equation_defects(hd, model, 2)) <= 2e-7
  @test maximum(native_hd_intrinsic_defects(hd, 2, model.d^2)) <= 2e-7
end

@testset "arbitrary native GKLS HD equals BF through order 4: full-rank sidebands" begin
  γ = 0.52
  model = NativeModel(
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
        -1 => sqrt(γ) * 0.31 * σx_native,
        0 => sqrt(γ) * 0.43 * σy_native,
        1 => sqrt(γ) * 0.57 * σz_native,
      ),
    ],
  )

  bf = native_bf_arbitrary(model, 4)
  hd = native_hd_arbitrary(model, 4)
  native_compare_arbitrary_states(bf, hd, 4)
  @test maximum(native_hd_equation_defects(hd, model, 4)) <= 5e-7
  @test maximum(native_hd_intrinsic_defects(hd, 4, model.d^2)) <= 5e-7

  ε = 2e-3
  @test norm(native_bf_finite(bf, model.d, ε) - native_hd_finite(hd, model.d, ε)) <= 1e-7
end

@testset "arbitrary native GKLS HD equals BF through order 4: Hamiltonian branch" begin
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

  bf = native_bf_arbitrary(model, 4)
  hd = native_hd_arbitrary(model, 4)
  native_compare_arbitrary_states(bf, hd, 4)

  @test isempty(hd.channels)
  @test all(norm(Sn) <= 1e-7 for Sn in hd.S)
  @test maximum(native_hd_equation_defects(hd, model, 4)) <= 5e-7
  @test maximum(native_hd_intrinsic_defects(hd, 4, model.d^2)) <= 5e-7
end

@testset "arbitrary native GKLS HD is prefix stable" begin
  γ = 0.48
  model = NativeModel(
    2,
    NativeFS(0 => 0.19 * σz_native, 1 => 0.27 * σx_native, -1 => 0.27 * σx_native),
    [NativeFS(0 => sqrt(γ) * σm_native)],
  )
  order2 = native_hd_arbitrary(model, 2)
  order3 = native_hd_arbitrary(model, 3)

  @test all(isapprox(order2.E[n], order3.E[n]; atol=2e-7, rtol=2e-7) for n in 1:3)
  @test all(isapprox(order2.S[n], order3.S[n]; atol=2e-7, rtol=2e-7) for n in 1:2)
end

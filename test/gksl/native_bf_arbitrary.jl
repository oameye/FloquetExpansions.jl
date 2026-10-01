include("native_bf_prototype.jl")

# Arbitrary-order native Bloch--Feshbach recurrence, run through the package driver
# FloquetExpansions.native_recurrence on the dense Liouville representation.
#
# The low-order oracle in native_bf_prototype.jl remains the independent frozen reference through
# generator order 2.

const NativeChannelSeries = FloquetExpansions.GradedChannel{ComplexF64}

struct NativeBFState
  E::Vector{NativeCM}
  S::Vector{NativeCM}
  bf_slots::Vector{NativeCM}
  H::Vector{NativeCM}
  Y::Vector{NativeFS}
  channels::Vector{NativeChannelSeries}
  known_gram::Vector{NativeCM}
end

native_arb_fsconst(X) = NativeFS(0 => NativeCM(X))

function native_arb_fsscale(A::NativeFS, α)
  return NativeFS(k => α * value for (k, value) in A)
end

function native_arb_series_zero(N)
  return [NativeFS() for _ in 0:N]
end

function native_arb_series_identity(n, N)
  out = native_arb_series_zero(N)
  out[1] = native_arb_fsconst(native_id(n))
  return out
end

function native_arb_series_add(A, B, N; α=1.0, β=1.0)
  return [native_fsadd(A[k], B[k], α, β) for k in 1:(N + 1)]
end

function native_arb_series_scale(A, N, α)
  return [native_arb_fsscale(A[k], α) for k in 1:(N + 1)]
end

function native_arb_series_mul(A, B, N)
  out = native_arb_series_zero(N)
  for n in 0:N
    coefficient = NativeFS()
    for j in 0:n
      j + 1 <= length(A) || continue
      n - j + 1 <= length(B) || continue
      coefficient = native_fsadd(coefficient, native_fsmul(A[j + 1], B[n - j + 1]))
    end
    out[n + 1] = coefficient
  end
  return out
end

function native_arb_series_exp(A, nsuper, N)
  result = native_arb_series_identity(nsuper, N)
  power = native_arb_series_identity(nsuper, N)
  for p in 1:N
    power = native_arb_series_mul(power, A, N)
    result = native_arb_series_add(
      result, native_arb_series_scale(power, N, inv(factorial(p))), N
    )
  end
  return result
end

function native_arb_series_log(K, nsuper, N)
  identity = native_arb_series_identity(nsuper, N)
  A = native_arb_series_add(K, identity, N; α=1.0, β=-1.0)
  result = native_arb_series_zero(N)
  power = A
  for p in 1:N
    coefficient = (-1.0)^(p + 1) / p
    result = native_arb_series_add(
      result, native_arb_series_scale(power, N, coefficient), N
    )
    p == N || (power = native_arb_series_mul(power, A, N))
  end
  return result
end

function native_arb_static_exp(S, nsuper, N; sign=-1.0)
  A = native_arb_series_zero(N)
  for n in 1:min(N, length(S))
    A[n + 1] = native_arb_fsconst(sign * S[n])
  end
  return native_arb_series_exp(A, nsuper, N)
end

function native_bf_intrinsic_offset(Y, Sprev, order, nsuper)
  @assert length(Y) >= order + 1
  Ytrunc = [copy(Y[k]) for k in 1:(order + 1)]
  # The yet-undetermined BF static slot at this order is set to zero.
  delete!(Ytrunc[order + 1], 0)
  minusS = native_arb_static_exp(Sprev, nsuper, order; sign=-1.0)
  Kosc = native_arb_series_mul(Ytrunc, minusS, order)
  logKosc = native_arb_series_log(Kosc, nsuper, order)
  return native_fsavg(logKosc[order + 1], nsuper)
end

function native_active_channels(channels, order, dim)
  return FloquetExpansions.active_channels(channels, order, dim)
end

native_known_gram(channels, order, dim) = FloquetExpansions.known_gram(channels, order, dim)

function native_store_active_corrections!(channels, indices, correction, order)
  return FloquetExpansions.store_corrections!(channels, indices, correction, order)
end

function native_gram_coefficient(channels, order, dim)
  return FloquetExpansions.gram_coefficient(channels, order, dim)
end

function native_set_static(A::NativeFS, value)
  result = copy(A)
  result[0] = value
  return result
end

function native_bf_static_residual(L, Y, E, order, Yn)
  residual = native_fsmul(L, Yn)
  for j in 1:order
    Yj = j == order ? Yn : Y[j + 1]
    residual = native_fsadd(
      residual, native_fs_right_static(Yj, E[order - j + 1]), 1.0, -1.0
    )
  end
  return residual
end

function native_bf_arbitrary(model::NativeModel, N; tol=1e-8, inverse=NativeNoInverse())
  result = FloquetExpansions.native_recurrence(
    FloquetExpansions.BlochFeshbach(),
    native_dense_representation(model.d),
    inverse,
    native_liouvillian_harmonics(model),
    native_sideband_columns(model),
    N,
    tol,
  )
  return NativeBFState(
    result.E,
    result.S,
    result.slots,
    result.H,
    result.kick,
    result.channels,
    result.known_gram,
  )
end

function native_bf_finite(state::NativeBFState, d, ε)
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

function native_bf_truncated(state::NativeBFState, ε)
  result = zero(state.E[1])
  for n in eachindex(state.E)
    result += ε^(n - 1) * state.E[n]
  end
  return result
end

function native_channel_prefix_equal(a, b, order; atol=1e-8)
  achannels = [channel for channel in a.channels if channel.onset <= order]
  bchannels = [channel for channel in b.channels if channel.onset <= order]
  length(achannels) == length(bchannels) || return false
  for (ca, cb) in zip(achannels, bchannels)
    ca.onset == cb.onset || return false
    maxcoefficient = order - ca.onset
    for k in 0:maxcoefficient
      k + 1 <= length(ca.coefficients) || return false
      k + 1 <= length(cb.coefficients) || return false
      isapprox(ca.coefficients[k + 1], cb.coefficients[k + 1]; atol, rtol=atol) ||
        return false
    end
  end
  return true
end

@testset "arbitrary native GKLS BF reproduces frozen orders 0-2" begin
  γ = 0.63
  model = NativeModel(
    2,
    NativeFS(0 => 0.29 * σz_native, 1 => 0.38 * σx_native, -1 => 0.38 * σx_native),
    [NativeFS(0 => sqrt(γ) * σm_native)],
  )
  oracle = native_bf_order02(model)
  result = native_bf_arbitrary(model, 2)

  @test all(isapprox(result.E[n], oracle.E[n]; atol=1e-8, rtol=1e-8) for n in 1:3)
  @test all(isapprox(result.S[n], oracle.S[n]; atol=1e-8, rtol=1e-8) for n in 1:2)
  @test all(
    isapprox(result.bf_slots[n], oracle.bf_slots[n]; atol=1e-8, rtol=1e-8) for n in 1:2
  )
  @test all(isapprox(result.H[n], oracle.H[n]; atol=1e-8, rtol=1e-8) for n in 1:3)
end

@testset "arbitrary native GKLS BF is prefix stable through order 4" begin
  γ = 0.52
  sideband_jump = NativeFS(
    -1 => sqrt(γ) * 0.31 * σx_native,
    0 => sqrt(γ) * 0.43 * σy_native,
    1 => sqrt(γ) * 0.57 * σz_native,
  )
  model = NativeModel(
    2,
    NativeFS(
      0 => 0.21 * σz_native,
      1 => 0.24 * σx_native,
      -1 => 0.24 * σx_native,
      2 => 0.07 * σy_native,
      -2 => 0.07 * σy_native,
    ),
    [sideband_jump],
  )

  order2 = native_bf_arbitrary(model, 2)
  order4 = native_bf_arbitrary(model, 4)

  @test all(isapprox(order2.E[n], order4.E[n]; atol=1e-8, rtol=1e-8) for n in 1:3)
  @test all(isapprox(order2.S[n], order4.S[n]; atol=1e-8, rtol=1e-8) for n in 1:2)
  @test native_channel_prefix_equal(order2, order4, 2)

  dim = model.d^2 - 1
  for order in 0:4
    gram = native_gram_coefficient(order4.channels, order, dim)
    @test isapprox(
      native_kossakowski(order4.E[order + 1], model.d), gram; atol=1e-7, rtol=1e-7
    )
  end

  ε = 2e-3
  finite = native_bf_finite(order4, model.d, ε)
  truncated = native_bf_truncated(order4, ε)
  @test norm(finite - truncated) <= 1e7 * ε^5
  @test norm(finite' * vec(native_id(model.d))) <= 1e-8
  cfinite = native_hermitian(native_kossakowski(finite, model.d))
  @test minimum(eigvals(Hermitian(cfinite))) >= -1e-9
end

@testset "a newborn channel joins the active flag at the next order" begin
  dim = 3
  leading = native_tlcoef(0.37 * σx_native, 2)
  channels = [NativeChannelSeries(1, [leading])]

  indices, active = native_active_channels(channels, 2, dim)
  @test indices == [1]
  @test active == reshape(leading, :, 1)
  @test iszero(native_known_gram(channels, 2, dim))

  correction = reshape(0.19 * native_tlcoef(σz_native, 2), :, 1)
  native_store_active_corrections!(channels, indices, correction, 2)
  expected2 = leading * correction[:, 1]' + correction[:, 1] * leading'
  @test isapprox(native_gram_coefficient(channels, 2, dim), expected2; atol=1e-12)

  expected3 = correction[:, 1] * correction[:, 1]'
  @test isapprox(native_known_gram(channels, 3, dim), expected3; atol=1e-12)
end

@testset "arbitrary native GKLS BF preserves the Hamiltonian branch through order 4" begin
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
  result = native_bf_arbitrary(model, 4)

  @test isempty(result.channels)
  @test all(norm(Sn) <= 1e-7 for Sn in result.S)
  @test all(norm(native_kossakowski(En, 2)) <= 1e-7 for En in result.E)
  @test all(
    norm(result.E[n] - native_lindblad(result.H[n], NativeCM[])) <= 1e-7 for
    n in eachindex(result.E)
  )
end

include("native_bf_prototype.jl")

# Research-only arbitrary-order Bloch--Feshbach recurrence for #350/#359.
#
# The low-order oracle in native_bf_prototype.jl remains the independent frozen reference through
# generator order 2. This file generalizes exactly the same semantics to arbitrary retained order.

mutable struct NativeChannelSeries
  onset::Int
  coefficients::Vector{Vector{ComplexF64}}
end

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
    result = native_arb_series_add(result, native_arb_series_scale(power, N, inv(factorial(p))), N)
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
    result = native_arb_series_add(result, native_arb_series_scale(power, N, coefficient), N)
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

function native_initial_channels(model::NativeModel)
  B0 = native_sideband_columns(model)
  return [NativeChannelSeries(0, [copy(B0[:, j])]) for j in axes(B0, 2)]
end

function native_active_channels(channels, order, dim)
  indices = Int[]
  columns = Vector{Vector{ComplexF64}}()
  for (index, channel) in pairs(channels)
    channel.onset < order || continue
    push!(indices, index)
    push!(columns, channel.coefficients[1])
  end
  active = isempty(columns) ? zeros(ComplexF64, dim, 0) : hcat(columns...)
  return indices, active
end

function native_known_gram(channels, order, dim)
  known = zeros(ComplexF64, dim, dim)
  for channel in channels
    channel.onset < order || continue
    q = order - channel.onset
    for k in 1:(q - 1)
      left = k + 1
      right = q - k + 1
      left <= length(channel.coefficients) || error("missing lower left amplitude coefficient")
      right <= length(channel.coefficients) || error("missing lower right amplitude coefficient")
      known += channel.coefficients[left] * channel.coefficients[right]'
    end
  end
  return native_hermitian(known)
end

function native_store_active_corrections!(channels, indices, correction, order)
  for (column, index) in pairs(indices)
    channel = channels[index]
    q = order - channel.onset
    length(channel.coefficients) == q || error("graded channel coefficient sequence is not prefix complete")
    push!(channel.coefficients, copy(correction[:, column]))
  end
  return channels
end

function native_store_births!(channels, newborn, order)
  for column in axes(newborn, 2)
    push!(channels, NativeChannelSeries(order, [copy(newborn[:, column])]))
  end
  return channels
end

function native_gram_coefficient(channels, order, dim)
  coefficient = zeros(ComplexF64, dim, dim)
  for channel in channels
    q = order - channel.onset
    q < 0 && continue
    for k in 0:q
      left = k + 1
      right = q - k + 1
      left <= length(channel.coefficients) || continue
      right <= length(channel.coefficients) || continue
      coefficient += channel.coefficients[left] * channel.coefficients[right]'
    end
  end
  return native_hermitian(coefficient)
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
      residual,
      native_fs_right_static(Yj, E[order - j + 1]),
      1.0,
      -1.0,
    )
  end
  return residual
end

function native_bf_arbitrary(model::NativeModel, N; tol=1e-8)
  N >= 0 || throw(ArgumentError("retained order must be nonnegative"))
  d = model.d
  nsuper = d^2
  dim = d^2 - 1
  L = native_liouvillian_harmonics(model)
  L0 = native_fsavg(L, nsuper)

  channels = native_initial_channels(model)
  E = NativeCM[L0]
  H = NativeCM[native_hamiltonian_part(L0, d)]
  S = NativeCM[]
  bf_slots = NativeCM[]
  known_gram = NativeCM[]

  c0 = native_gram_coefficient(channels, 0, dim)
  norm(native_kossakowski(L0, d) - c0) <= tol * max(1.0, norm(c0)) ||
    error("order-zero native channels do not reconstruct the averaged Kossakowski tensor")

  Y = NativeFS[native_arb_fsconst(native_id(nsuper))]
  N == 0 && return NativeBFState(E, S, bf_slots, H, Y, channels, known_gram)
  push!(Y, native_fsint(L))

  for order in 1:N
    offset = native_bf_intrinsic_offset(Y, S, order, nsuper)
    base_slot = -offset
    Ynbase = native_set_static(Y[order + 1], base_slot)
    residual_base = native_bf_static_residual(L, Y, E, order, Ynbase)
    Vhat = native_fsavg(residual_base, nsuper)

    known = native_known_gram(channels, order, dim)
    active_indices, active = native_active_channels(channels, order, dim)
    step = native_static_step(L0, Vhat, known, active, d; tol)

    push!(S, step.S)
    push!(bf_slots, base_slot + step.S)
    push!(E, step.E)
    push!(H, step.H)
    push!(known_gram, known)

    Y[order + 1] = native_set_static(Y[order + 1], bf_slots[end])
    native_store_active_corrections!(channels, active_indices, step.correction, order)
    native_store_births!(channels, step.newborn, order)

    reconstructed = native_gram_coefficient(channels, order, dim)
    scale = max(1.0, norm(reconstructed))
    norm(native_kossakowski(E[end], d) - reconstructed) <= 20 * tol * scale ||
      error("graded channel state failed coefficient reconstruction")

    if order < N
      forcing = native_bf_static_residual(L, Y, E, order, Y[order + 1])
      push!(Y, native_fsint(forcing))
    end
  end

  return NativeBFState(E, S, bf_slots, H, Y, channels, known_gram)
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
      isapprox(ca.coefficients[k + 1], cb.coefficients[k + 1]; atol, rtol=atol) || return false
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
    @test isapprox(native_kossakowski(order4.E[order + 1], model.d), gram; atol=1e-7, rtol=1e-7)
  end

  ε = 2e-3
  finite = native_bf_finite(order4, model.d, ε)
  truncated = native_bf_truncated(order4, ε)
  @test norm(finite - truncated) <= 1e7 * ε^5
  @test norm(finite' * vec(native_id(model.d))) <= 1e-8
  cfinite = native_hermitian(native_kossakowski(finite, model.d))
  @test minimum(eigvals(Hermitian(cfinite))) >= -1e-9
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
    norm(result.E[n] - native_lindblad(result.H[n], NativeCM[])) <= 1e-7 for n in eachindex(result.E)
  )
end

using Test
using LinearAlgebra

# Research-only dense-matrix oracle for #350/#352.
#
# This intentionally does not touch the public API. It implements the direct Bloch--Feshbach
# recurrence only through generator order 2 and solves the retained m=0 coefficient in GKLS
# coordinates before advancing the oscillatory recurrence.

const NativeCM = Matrix{ComplexF64}
const NativeFS = Dict{Int,NativeCM}

struct NativeModel
  d::Int
  H::NativeFS
  jumps::Vector{NativeFS}
end

struct NativeOrder02
  E::Vector{NativeCM}
  S::Vector{NativeCM}
  bf_slots::Vector{NativeCM}
  H::Vector{NativeCM}
  birth0::Vector{NativeCM}
  birth1::Vector{NativeCM}
  birth2::Vector{NativeCM}
  Y1osc::NativeFS
  Y2osc::NativeFS
  K2known::NativeCM
end

native_id(d) = Matrix{ComplexF64}(I, d, d)
native_lmul(A) = kron(native_id(size(A, 1)), A)
native_rmul(B) = kron(transpose(B), native_id(size(B, 1)))
native_sand(X, Y) = kron(conj(Y), X)

function native_lindblad(H, jumps)
  result = -im * (native_lmul(H) - native_rmul(H))
  for R in jumps
    norm2 = R' * R
    result +=
      native_sand(R, R) - 0.5 * native_lmul(norm2) - 0.5 * native_rmul(norm2)
  end
  return result
end

function native_fsadd(A::NativeFS, B::NativeFS, α=1.0, β=1.0)
  C = NativeFS()
  for (k, value) in A
    C[k] = α * value
  end
  for (k, value) in B
    C[k] = get(C, k, zero(value)) + β * value
  end
  return C
end

function native_fsmul(A::NativeFS, B::NativeFS)
  C = NativeFS()
  for (ka, a) in A, (kb, b) in B
    k = ka + kb
    C[k] = get(C, k, zero(a * b)) + a * b
  end
  return C
end

native_fsavg(A::NativeFS, n) = get(A, 0, zeros(ComplexF64, n, n))
native_fsint(A::NativeFS) =
  NativeFS(k => (im / k) * value for (k, value) in A if k != 0)

function native_fs_right_static(A::NativeFS, B::NativeCM)
  return NativeFS(k => value * B for (k, value) in A)
end

function native_liouvillian_harmonics(model::NativeModel)
  d = model.d
  out = NativeFS()
  function add!(m, value)
    out[m] = get(out, m, zeros(ComplexF64, d^2, d^2)) + value
  end
  for (m, Hm) in model.H
    add!(m, -im * (native_lmul(Hm) - native_rmul(Hm)))
  end
  for jump in model.jumps, (p, Rp) in jump, (q, Rq) in jump
    normpq = Rq' * Rp
    add!(
      p - q,
      native_sand(Rp, Rq) - 0.5 * native_lmul(normpq) - 0.5 * native_rmul(normpq),
    )
  end
  return out
end

function native_traceless_basis(d)
  basis = NativeCM[]
  for j in 1:d, k in (j + 1):d
    X = zeros(ComplexF64, d, d)
    X[j, k] = 1
    X[k, j] = 1
    push!(basis, X / sqrt(2))

    Y = zeros(ComplexF64, d, d)
    Y[j, k] = -im
    Y[k, j] = im
    push!(basis, Y / sqrt(2))
  end
  for l in 1:(d - 1)
    Z = zeros(ComplexF64, d, d)
    for j in 1:l
      Z[j, j] = 1
    end
    Z[l + 1, l + 1] = -l
    push!(basis, Z / sqrt(l * (l + 1)))
  end
  return basis
end

function native_full_basis(d)
  return vcat([native_id(d) / sqrt(d)], native_traceless_basis(d))
end

function native_chimat(L, d)
  basis = native_full_basis(d)
  n = d^2
  χ = zeros(ComplexF64, n, n)
  for i in 1:n, j in 1:n
    χ[i, j] = dot(vec(native_sand(basis[i], basis[j])), vec(L))
  end
  return χ
end

native_kossakowski(L, d) = native_chimat(L, d)[2:end, 2:end]
native_hermitian(X) = (X + X') / 2

function native_hamiltonian_part(L, d)
  χ = native_chimat(L, d)
  basis = native_traceless_basis(d)
  H = zeros(ComplexF64, d, d)
  for i in eachindex(basis)
    H += χ[i + 1, 1] / sqrt(d) * basis[i]
  end
  return native_hermitian((im / 2) * (H - H'))
end

function native_from_Hc(H, c, d)
  basis = native_traceless_basis(d)
  result = -im * (native_lmul(H) - native_rmul(H))
  for i in eachindex(basis), j in eachindex(basis)
    abs(c[i, j]) <= 1e-14 && continue
    Fji = basis[j]' * basis[i]
    result += c[i, j] * (
      native_sand(basis[i], basis[j]) -
      0.5 * native_lmul(Fji) -
      0.5 * native_rmul(Fji)
    )
  end
  return result
end

function native_tlcoef(X, d)
  return ComplexF64[tr(F' * X) for F in native_traceless_basis(d)]
end

function native_operator(v, d)
  basis = native_traceless_basis(d)
  result = zeros(ComplexF64, d, d)
  for i in eachindex(basis)
    result += v[i] * basis[i]
  end
  return result
end

function native_hvec(X)
  n = size(X, 1)
  values = Float64[]
  for i in 1:n
    push!(values, real(X[i, i]))
  end
  for i in 1:n, j in (i + 1):n
    push!(values, sqrt(2) * real(X[i, j]))
    push!(values, sqrt(2) * imag(X[i, j]))
  end
  return values
end

function native_hmat(v, n)
  X = zeros(ComplexF64, n, n)
  k = 1
  for i in 1:n
    X[i, i] = v[k]
    k += 1
  end
  for i in 1:n, j in (i + 1):n
    X[i, j] = (v[k] + im * v[k + 1]) / sqrt(2)
    X[j, i] = conj(X[i, j])
    k += 2
  end
  return X
end

const NATIVE_GAUGE_CACHE = Dict{Int,Vector{NativeCM}}()

function native_gauge_algebra(d)
  return get!(NATIVE_GAUGE_CACHE, d) do
    basis = native_full_basis(d)
    n = d^2

    function superop(χ)
      result = zeros(ComplexF64, n, n)
      for i in 1:n, j in 1:n
        result += χ[i, j] * native_sand(basis[i], basis[j])
      end
      return result
    end

    function tp_operator(χ)
      result = zeros(ComplexF64, d, d)
      for i in 1:n, j in 1:n
        result += χ[i, j] * basis[j]' * basis[i]
      end
      return native_hermitian(result)
    end

    columns = Vector{Vector{Float64}}()
    for k in 1:(n^2)
      χ = native_hmat(Float64.(1:(n^2) .== k), n)
      push!(columns, native_hvec(tp_operator(χ)))
    end
    constraints = hcat(columns...)
    kernel = nullspace(constraints)
    return [superop(native_hmat(kernel[:, j], n)) for j in axes(kernel, 2)]
  end
end

function native_sideband_columns(model::NativeModel; tol=1e-12)
  columns = Vector{Vector{ComplexF64}}()
  for jump in model.jumps
    for (_, R) in sort!(collect(jump); by=first)
      coefficient = native_tlcoef(R, model.d)
      norm(coefficient) > tol && push!(columns, coefficient)
    end
  end
  n = model.d^2 - 1
  if isempty(columns)
    return zeros(ComplexF64, n, 0)
  end
  return hcat(columns...)
end

function native_active_split(A; tol=1e-10)
  n = size(A, 1)
  if size(A, 2) == 0
    return zeros(ComplexF64, n, 0), native_id(n)
  end
  decomposition = svd(A; full=true)
  scale = max(1.0, maximum(decomposition.S))
  rankA = count(>(tol * scale), decomposition.S)
  U = decomposition.U[:, 1:rankA]
  P = decomposition.U[:, (rankA + 1):n]
  return U, P
end

function native_phi_matrix(L0, d, P, gauge_basis)
  darkdim = size(P, 2)
  if darkdim == 0
    return zeros(Float64, 0, length(gauge_basis))
  end
  columns = Vector{Vector{Float64}}()
  for S in gauge_basis
    commutator = L0 * S - S * L0
    dark = native_hermitian(P' * native_kossakowski(commutator, d) * P)
    push!(columns, native_hvec(dark))
  end
  if isempty(columns)
    return zeros(Float64, darkdim^2, 0)
  end
  return hcat(columns...)
end

function native_dark_solve(L0, Vhat, known, active, d; tol=1e-9)
  _, P = native_active_split(active; tol)
  darkdim = size(P, 2)
  gauge_basis = native_gauge_algebra(d)
  if darkdim == 0
    return (
      S=zeros(ComplexF64, d^2, d^2),
      newborn=zeros(ComplexF64, d^2 - 1, 0),
      dark_residual=zeros(ComplexF64, 0, 0),
    )
  end

  delta = native_hermitian(P' * (native_kossakowski(Vhat, d) - known) * P)
  Φ = native_phi_matrix(L0, d, P, gauge_basis)
  if isempty(gauge_basis)
    coordinates = Float64[]
  else
    coordinates = -pinv(Φ; rtol=1e-10) * native_hvec(delta)
  end
  S = zeros(ComplexF64, d^2, d^2)
  for i in eachindex(coordinates)
    S += coordinates[i] * gauge_basis[i]
  end

  residual = native_hermitian(
    delta + P' * native_kossakowski(L0 * S - S * L0, d) * P
  )
  scale = max(1.0, norm(delta))
  eig = eigen(Hermitian(residual))
  minimum(eig.values) >= -tol * scale || error(
    "order-0/1/2 prototype found a non-PSD dark cokernel; full cone solve is not implemented",
  )

  keep = findall(>(tol * scale), eig.values)
  if isempty(keep)
    newborn = zeros(ComplexF64, d^2 - 1, 0)
  else
    newborn = P * eig.vectors[:, keep] * Diagonal(sqrt.(eig.values[keep]))
  end
  return (; S, newborn, dark_residual=residual)
end

function native_tangent_lift(active, target; tol=1e-9)
  n, r = size(active)
  if r == 0
    norm(target) <= tol || error("nonzero tangent target with no active channels")
    return zeros(ComplexF64, n, 0)
  end

  columns = Vector{Vector{Float64}}()
  for j in 1:r, i in 1:n
    for phase in (1.0 + 0im, 0.0 + 1im)
      D = zeros(ComplexF64, n, r)
      D[i, j] = phase
      push!(columns, native_hvec(native_hermitian(active * D' + D * active')))
    end
  end
  map = hcat(columns...)
  coordinates = pinv(map; rtol=1e-10) * native_hvec(native_hermitian(target))
  D = zeros(ComplexF64, n, r)
  k = 1
  for j in 1:r, i in 1:n
    D[i, j] = coordinates[k] + im * coordinates[k + 1]
    k += 2
  end
  residual = target - active * D' - D * active'
  norm(residual) <= tol * max(1.0, norm(target)) || error("tangent Gram lift failed")
  return D
end

function native_static_step(L0, Vhat, known, active, d; tol=1e-8)
  dark = native_dark_solve(L0, Vhat, known, active, d; tol)
  E = Vhat + L0 * dark.S - dark.S * L0
  Pfull = dark.newborn * dark.newborn'
  target = native_hermitian(native_kossakowski(E, d) - known - Pfull)
  correction = native_tangent_lift(active, target; tol)
  reconstructed = known + active * correction' + correction * active' + Pfull
  scale = max(1.0, norm(reconstructed))
  norm(native_kossakowski(E, d) - reconstructed) <= tol * scale ||
    error("native static step failed to reconstruct its Kossakowski coefficient")
  H = native_hamiltonian_part(E, d)
  norm(E - native_from_Hc(H, reconstructed, d)) <= 20 * tol * max(1.0, norm(E)) ||
    error("native static step failed Hamiltonian/Kossakowski reconstruction")
  return (; E, S=dark.S, H, correction, newborn=dark.newborn, C=reconstructed)
end

function native_bf_order02(model::NativeModel; tol=1e-8)
  d = model.d
  nsuper = d^2
  L = native_liouvillian_harmonics(model)
  L0 = native_fsavg(L, nsuper)

  # Order 0.
  E0 = L0
  Y1osc = native_fsint(L)
  B0 = native_sideband_columns(model)
  c0 = B0 * B0'
  norm(native_kossakowski(E0, d) - c0) <= tol * max(1.0, norm(c0)) ||
    error("order-zero sideband columns do not reconstruct the averaged Kossakowski tensor")
  H0 = native_hamiltonian_part(E0, d)

  # Order 1: solve the static GKLS slot before computing Y2^osc.
  V1 = native_fsavg(native_fsmul(L, Y1osc), nsuper)
  known1 = zeros(ComplexF64, d^2 - 1, d^2 - 1)
  step1 = native_static_step(L0, V1, known1, B0, d; tol)
  S1 = step1.S
  E1 = step1.E
  B1 = step1.correction
  U0 = step1.newborn

  Y1 = copy(Y1osc)
  Y1[0] = S1
  forcing1 = native_fsadd(
    native_fsmul(L, Y1), native_fs_right_static(Y1, E0), 1.0, -1.0
  )
  Y2osc = native_fsint(forcing1)

  # Order 2 in intrinsic total-kick coordinates.
  V2 = native_fsavg(native_fsmul(L, Y2osc), nsuper) - S1 * E1
  Y1square_avg = native_fsavg(native_fsmul(Y1osc, Y1osc), nsuper)
  a2bf = -0.5 * (S1 * S1 + Y1square_avg)
  Vhat2 = V2 - (L0 * a2bf - a2bf * L0)

  K2known = B1 * B1'
  active2 = hcat(B0, U0)
  step2 = native_static_step(L0, Vhat2, K2known, active2, d; tol)
  S2 = step2.S
  E2 = step2.E
  correction2 = step2.correction
  B2 = correction2[:, 1:size(B0, 2)]
  U1 = correction2[:, (size(B0, 2) + 1):end]
  V0 = step2.newborn
  s2bf = S2 - a2bf

  # Direct BF coordinate identity at order two.
  E2bf = V2 + L0 * s2bf - s2bf * L0
  norm(E2 - E2bf) <= tol * max(1.0, norm(E2)) ||
    error("intrinsic/BF order-two slots disagree")

  # birth0[j] = [R_j,0, R_j,1, R_j,2], birth1[j] = [U_j,0, U_j,1], and
  # birth2[j] = [V_j,0]. The vectors are flattened to keep this fixture local to the tests.
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

  return NativeOrder02(
    [E0, E1, E2],
    [S1, S2],
    [S1, s2bf],
    [H0, step1.H, step2.H],
    birth0,
    birth1,
    birth2,
    Y1osc,
    Y2osc,
    K2known,
  )
end

function native_effective_gkls(result::NativeOrder02, ε)
  H = result.H[1] + ε * result.H[2] + ε^2 * result.H[3]
  jumps = NativeCM[]
  for k in 1:3:length(result.birth0)
    push!(
      jumps,
      result.birth0[k] + ε * result.birth0[k + 1] + ε^2 * result.birth0[k + 2],
    )
  end
  for k in 1:2:length(result.birth1)
    push!(jumps, sqrt(ε) * (result.birth1[k] + ε * result.birth1[k + 1]))
  end
  for R in result.birth2
    push!(jumps, ε * R)
  end
  return native_lindblad(H, jumps)
end

σx_native = ComplexF64[0 1; 1 0]
σy_native = ComplexF64[0 -im; im 0]
σz_native = ComplexF64[1 0; 0 -1]
σm_native = ComplexF64[0 0; 1 0]

@testset "native GKLS BF research prototype: Hamiltonian reduction through order 2" begin
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
  result = native_bf_order02(model)

  @test norm(result.S[1]) <= 1e-9
  @test norm(result.S[2]) <= 1e-9
  @test isempty(result.birth0)
  @test isempty(result.birth1)
  @test isempty(result.birth2)
  @test all(norm(native_kossakowski(result.E[n], 2)) <= 1e-9 for n in 1:3)
  @test all(
    norm(result.E[n] - native_lindblad(result.H[n], NativeCM[])) <= 1e-8 for n in 1:3
  )
end

@testset "native GKLS BF research prototype: static-loss dark solve" begin
  γ = 0.63
  model = NativeModel(
    2,
    NativeFS(
      0 => 0.29 * σz_native,
      1 => 0.38 * σx_native,
      -1 => 0.38 * σx_native,
    ),
    [NativeFS(0 => sqrt(γ) * σm_native)],
  )
  result = native_bf_order02(model)

  @test length(result.birth0) == 3
  @test norm(result.E[1] - native_lindblad(result.H[1], [result.birth0[1]])) <= 1e-8
  @test all(isfinite, result.S[1])
  @test all(isfinite, result.S[2])

  ε = 1e-3
  finite = native_effective_gkls(result, ε)
  truncated = result.E[1] + ε * result.E[2] + ε^2 * result.E[3]
  @test norm(finite - truncated) <= 1e3 * ε^3
  cfinite = native_hermitian(native_kossakowski(finite, 2))
  @test minimum(eigvals(Hermitian(cfinite))) >= -1e-10
end

@testset "native GKLS BF research prototype: sidebands are native columns" begin
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
    ),
    [sideband_jump],
  )
  result = native_bf_order02(model)

  # Three physical jump sidebands span the full qubit traceless space at birth order zero.
  @test length(result.birth0) == 9
  @test isempty(result.birth1)
  @test isempty(result.birth2)

  B1 = hcat(
    native_tlcoef(result.birth0[2], 2),
    native_tlcoef(result.birth0[5], 2),
    native_tlcoef(result.birth0[8], 2),
  )
  @test norm(result.K2known - B1 * B1') <= 1e-9

  # This is the low-order reason no external transported-channel allowance Γ is needed:
  # the square of the first sideband-amplitude correction is ordinary K_2^< state.
  @test norm(result.K2known) >= 0
end

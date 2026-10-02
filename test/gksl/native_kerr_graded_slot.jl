include("native_kerr_fock.jl")
include("native_kerr_localized.jl")

# Research-only resolvent section of the native static slot for U(1)-symmetric averaged
# generators.
#
# ad_Lbar preserves the superoperator charge Q = (out_L - in_L) - (out_R - in_R). On the charged
# sectors Q != 0 of driven Kerr with loss every Bohr multiplier is nonzero on the physical lattice,
# so ad_Lbar is invertible there. The slot then chooses S_n with
#
#   [Lbar, S_n] = -Pi_{Q != 0} D[pi Delta_n pi],
#
# which cancels exactly the charged dark-dark target and touches nothing else. Unlike the
# minimum-norm section this needs no metric on the gauge algebra, and it is precisely what the
# localized ladder (D + J)^{-1} of #376/#377 computes. The neutral sector Q = 0 keeps the dense
# affine PSD section, where the newborn Gram form P_n lives.

function native_superoperator_charge(d)
  Q = zeros(Int, d^2, d^2)
  for output in 1:(d ^ 2), input in 1:(d ^ 2)
    outleft, outright = (output - 1) % d, (output - 1) ÷ d
    inleft, inright = (input - 1) % d, (input - 1) ÷ d
    Q[output, input] = (outleft - inleft) - (outright - inright)
  end
  return Q
end

function native_superoperator_level(d)
  level = zeros(Int, d^2, d^2)
  for output in 1:(d ^ 2), input in 1:(d ^ 2)
    level[output, input] = (output - 1) ÷ d - (input - 1) ÷ d
  end
  return level
end

native_charged_part(X, Q) = X .* (Q .!= 0)
native_neutral_part(X, Q) = X .* (Q .== 0)

const NATIVE_NEUTRAL_GAUGE_CACHE = Dict{Int,Vector{NativeCM}}()

function native_neutral_gauge_algebra(d)
  return get!(NATIVE_NEUTRAL_GAUGE_CACHE, d) do
    Q = native_superoperator_charge(d)
    projected = [native_neutral_part(G, Q) for G in native_gauge_algebra(d)]
    # Real coefficients keep every basis element Hermiticity preserving.
    M = hcat([vcat(real(vec(G)), imag(vec(G))) for G in projected]...)
    F = svd(M)
    rank = count(>(1e-10 * max(1.0, F.S[1])), F.S)
    return [
      sum(F.V[j, i] * projected[j] for j in eachindex(projected)) / F.S[i] for i in 1:rank
    ]
  end
end

function native_charged_dark_target(Vhat, known, P, d)
  Q = native_superoperator_charge(d)
  size(P, 2) == 0 && return zeros(ComplexF64, d^2, d^2)
  delta = native_hermitian(P' * (native_kossakowski(Vhat, d) - known) * P)
  dissipator = native_from_Hc(zeros(ComplexF64, d, d), P * delta * P', d)
  return -native_charged_part(dissipator, Q)
end

function native_charged_dark_kossakowski(C, known, active, d)
  Q = native_superoperator_charge(d)
  P = FloquetExpansions.gram_active_frame(Matrix{ComplexF64}(active); rtol=1e-8).dark
  size(P, 2) == 0 && return zeros(ComplexF64, d^2, d^2)
  delta = native_hermitian(P' * (C - known) * P)
  dissipator = native_from_Hc(zeros(ComplexF64, d, d), P * delta * P', d)
  return -native_charged_part(dissipator, Q)
end

function native_checked_charged_solve(inverse, L0, Y, tol)
  S = native_charged_solve(inverse, Y)
  defect = norm(L0 * S - S * L0 - Y)
  defect <= tol * max(1.0, norm(Y)) ||
    error("charged homological inverse does not solve [Lbar,S]=Y: defect $defect")
  return S
end

# Generic finite-dimensional reference: dense LU of ad_Lbar restricted to all Q != 0 entries.
struct NativeDenseChargedInverse
  d::Int
  entries::Vector{Tuple{Int,Int}}
  factor::LU{ComplexF64,Matrix{ComplexF64},Vector{Int}}
end

function NativeDenseChargedInverse(L0, d)
  Q = native_superoperator_charge(d)
  entries = [
    (output, input) for input in 1:(d ^ 2), output in 1:(d ^ 2) if Q[output, input] != 0
  ]
  m = length(entries)
  A = zeros(ComplexF64, m, m)
  for (column, (output, input)) in pairs(entries), (row, (rowout, rowin)) in pairs(entries)
    value = zero(ComplexF64)
    rowin == input && (value += L0[rowout, output])
    rowout == output && (value -= L0[input, rowin])
    A[row, column] = value
  end
  return NativeDenseChargedInverse(d, entries, lu(A))
end

function native_charged_solve(inverse::NativeDenseChargedInverse, Y)
  x = inverse.factor \ ComplexF64[Y[output, input] for (output, input) in inverse.entries]
  S = zeros(ComplexF64, inverse.d^2, inverse.d^2)
  for (k, (output, input)) in pairs(inverse.entries)
    S[output, input] = x[k]
  end
  return S
end

# Algebra-specific inverse: the localized affine Kerr Bohr multipliers d_t(n_L, n_R) on the
# diagonal and the triangular recycling ladder s_t = (y_t - J s_{t+1}) / d_t, with no linear solve.
struct NativeKerrLadderInverse
  d::Int
  κ::Float64
  charge::Matrix{Int}
  level::Matrix{Int}
  bohr::Matrix{ComplexF64}
  recycling::NativeCM
end

function NativeKerrLadderInverse(d, Δ, χ, κ; tol=1e-10)
  charge = native_superoperator_charge(d)
  level = native_superoperator_level(d)
  bohr = zeros(ComplexF64, d^2, d^2)
  for output in 1:(d ^ 2), input in 1:(d ^ 2)
    charge[output, input] == 0 && continue
    inleft, inright = (input - 1) % d, (input - 1) ÷ d
    qleft = (output - 1) % d - inleft
    qright = (output - 1) ÷ d - inright
    multiplier = localized_kerr_bohr2(qleft, qright, Δ, χ, κ)
    bohr[output, input] = ComplexF64(localized_value(multiplier, inleft, inright))
    abs(bohr[output, input]) > tol ||
      error("physical Kerr resonance in charged static sector")
  end
  a = native_fock_annihilation(d)
  return NativeKerrLadderInverse(d, κ, charge, level, bohr, native_sand(a, a))
end

function native_charged_solve(inverse::NativeKerrLadderInverse, Y)
  S = zeros(ComplexF64, inverse.d^2, inverse.d^2)
  A = inverse.recycling
  for t in (inverse.d - 1):-1:(1 - inverse.d)
    rung = (inverse.level .== t) .& (inverse.charge .!= 0)
    # Recycling maps rung t + 1 onto rung t, so S already holds every rung that feeds this one.
    rhs = Y - inverse.κ * (A * S - S * A)
    S[rung] = rhs[rung] ./ inverse.bohr[rung]
  end
  return S
end

struct NativeChargedSection{T} <: FloquetExpansions.HomologicalInverse
  inverse::T
  L0::NativeCM
  tol::Float64
end

function FloquetExpansions.regular_gauge(
  section::NativeChargedSection,
  representation::FloquetExpansions.DenseLiouvilleRepresentation,
  L0,
  C,
  known,
  active,
  tol,
)
  d = representation.d
  norm(L0 - section.L0) <= section.tol * max(1.0, norm(L0)) ||
    error("charged section was built for a different averaged generator")
  Y = native_charged_dark_kossakowski(C, known, active, d)
  return native_checked_charged_solve(section.inverse, L0, Y, section.tol)
end

function FloquetExpansions.singular_gauge_directions(
  ::NativeChargedSection, representation::FloquetExpansions.DenseLiouvilleRepresentation
)
  return native_neutral_gauge_algebra(representation.d)
end

function native_kerr_model(d, Δ, χ, κ, F)
  a = native_fock_annihilation(d)
  n = a' * a
  x = a + a'
  jumps = iszero(κ) ? NativeFS[] : [NativeFS(0 => sqrt(κ) * a)]
  return NativeModel(d, NativeFS(0 => Δ * n + χ * n * n, 1 => F * x, -1 => F' * x), jumps)
end

function native_ladder_section(model, Δ, χ, κ)
  L0 = native_fsavg(native_liouvillian_harmonics(model), model.d^2)
  return NativeChargedSection(NativeKerrLadderInverse(model.d, Δ, χ, κ), L0, 1e-9)
end

function native_dense_section(model)
  L0 = native_fsavg(native_liouvillian_harmonics(model), model.d^2)
  return NativeChargedSection(NativeDenseChargedInverse(L0, model.d), L0, 1e-9)
end

function native_sorted_spectrum(X)
  return sort(eigvals(X); by=z -> (round(real(z); digits=6), imag(z)))
end

const NATIVE_KERR = (Δ=1 / 2, χ=3 / 10, κ=4 / 5, F=2 / 5 + im / 5)

@testset "Kerr ladder inverse equals the dense charged inverse" begin
  (; Δ, χ, κ) = NATIVE_KERR
  for d in (3, 4)
    model = native_kerr_model(d, Δ, χ, κ, NATIVE_KERR.F)
    L0 = native_fsavg(native_liouvillian_harmonics(model), d^2)
    Q = native_superoperator_charge(d)
    @test norm(native_charged_part(L0, Q)) == 0

    ladder = NativeKerrLadderInverse(d, Δ, χ, κ)
    dense = NativeDenseChargedInverse(L0, d)
    G = native_gauge_algebra(d)
    Y = native_charged_part(L0 * G[3] - G[3] * L0 + 0.7 * (L0 * G[end] - G[end] * L0), Q)
    Sladder = native_checked_charged_solve(ladder, L0, Y, 1e-10)
    Sdense = native_checked_charged_solve(dense, L0, Y, 1e-10)
    @test norm(Sladder - Sdense) <= 1e-10 * max(1.0, norm(Sdense))
    @test norm(native_neutral_part(Sladder, Q)) == 0
  end
end

@testset "a physical Kerr resonance is rejected by the ladder" begin
  # Delta + chi (n_L + n_R) = 0 at n_L + n_R = 1 makes d_t vanish in the Q = 2 sector.
  @test_throws ErrorException NativeKerrLadderInverse(3, -0.3, 0.3, 0.8)
end

@testset "graded native slot through order 4 on driven Kerr" begin
  (; Δ, χ, κ, F) = NATIVE_KERR
  d = 3
  model = native_kerr_model(d, Δ, χ, κ, F)
  ladder = native_ladder_section(model, Δ, χ, κ)
  dense = native_dense_section(model)

  bf = native_bf_arbitrary(model, 4; inverse=ladder)
  bfdense = native_bf_arbitrary(model, 4; inverse=dense)
  hd = native_hd_arbitrary(model, 4; inverse=ladder)
  bf2 = native_bf_arbitrary(model, 2; inverse=ladder)
  minimum_norm = native_bf_arbitrary(model, 4)

  # The algebra-specific ladder and the generic dense inverse define the same section.
  @test all(norm(bf.E[k] - bfdense.E[k]) <= 1e-8 * max(1.0, norm(bf.E[k])) for k in 1:5)
  @test all(norm(bf.S[k] - bfdense.S[k]) <= 1e-8 * max(1.0, norm(bf.S[k])) for k in 1:4)
  @test native_channel_prefix_equal(bf, bfdense, 4; atol=1e-8)

  # The section is intrinsic, so the independent HD recurrence reproduces it.
  native_compare_arbitrary_states(bf, hd, 4; tol=2e-6)
  @test maximum(native_hd_equation_defects(hd, model, 4)) <= 2e-6
  @test maximum(native_hd_intrinsic_defects(hd, 4, d^2)) <= 2e-6

  @test all(isapprox(bf2.E[k], bf.E[k]; atol=2e-7, rtol=2e-7) for k in 1:3)
  @test native_channel_prefix_equal(bf2, bf, 2; atol=2e-7)

  dim = d^2 - 1
  Q = native_superoperator_charge(d)
  charged_orders = Int[]
  for order in 1:4
    @test isapprox(
      native_kossakowski(bf.E[order + 1], d),
      native_gram_coefficient(bf.channels, order, dim);
      atol=2e-6,
      rtol=2e-6,
    )
    norm(native_charged_part(bf.S[order], Q)) > 1e-8 && push!(charged_orders, order)

    # Defining property: no charged dark-dark target survives in E_n.
    _, active = native_active_channels(bf.channels, order, dim)
    _, P = native_active_split(active; tol=1e-8)
    residual = native_charged_dark_target(bf.E[order + 1], bf.known_gram[order], P, d)
    @test norm(residual) <= 1e-8 * max(1.0, norm(bf.E[order + 1]))
  end
  @test charged_orders == [4]

  # Through order 3 the charged sector is empty, so the resolvent and minimum-norm sections
  # coincide. At order 4 they differ by a static gauge in ker Phi_4: same newborn Gram form, same
  # dark-dark coefficient, different representative.
  @test all(isapprox(bf.E[k], minimum_norm.E[k]; atol=1e-8, rtol=1e-8) for k in 1:4)
  @test norm(bf.S[4] - minimum_norm.S[4]) > 1e-3
  _, active4 = native_active_channels(bf.channels, 4, dim)
  _, P4 = native_active_split(active4; tol=1e-8)
  dark(E) = P4' * native_kossakowski(E, d) * P4
  @test norm(dark(bf.E[5]) - dark(minimum_norm.E[5])) <= 1e-8
  δS = bf.S[4] - minimum_norm.S[4]
  # Phi_n is charge graded, so the neutral components of both sections agree.
  @test norm(native_neutral_part(δS, Q)) <= 1e-10
  @test norm(bf.E[5] - minimum_norm.E[5] - (bf.E[1] * δS - δS * bf.E[1])) <= 1e-8

  # Both sections lie on the same static Floquet orbit, so their truncated generators are
  # similar up to O(ε^5).
  defects = map((2e-2, 1e-2)) do ε
    a = native_sorted_spectrum(native_bf_truncated(bf, ε))
    b = native_sorted_spectrum(native_bf_truncated(minimum_norm, ε))
    maximum(abs, a - b)
  end
  @test defects[2] <= 1e-9
  @test defects[2] <= defects[1] / 24

  ε = 5e-4
  finite = native_bf_finite(bf, d, ε)
  @test norm(finite - native_hd_finite(hd, d, ε)) <= 2e-6
  @test norm(finite' * vec(native_id(d))) <= 2e-8
  @test minimum(eigvals(Hermitian(native_hermitian(native_kossakowski(finite, d))))) >=
    -2e-8
end

@testset "graded native slot at a larger Fock cutoff" begin
  (; Δ, χ, κ, F) = NATIVE_KERR
  d = 4
  model = native_kerr_model(d, Δ, χ, κ, F)
  bf = native_bf_arbitrary(model, 4; inverse=native_ladder_section(model, Δ, χ, κ))
  bfdense = native_bf_arbitrary(model, 4; inverse=native_dense_section(model))
  @test all(norm(bf.E[k] - bfdense.E[k]) <= 1e-8 * max(1.0, norm(bf.E[k])) for k in 1:5)
  @test native_channel_prefix_equal(bf, bfdense, 4; atol=1e-8)
  @test norm(native_charged_part(bf.S[4], native_superoperator_charge(d))) > 1e-3
  for order in 0:4
    @test isapprox(
      native_kossakowski(bf.E[order + 1], d),
      native_gram_coefficient(bf.channels, order, d^2 - 1);
      atol=2e-6,
      rtol=2e-6,
    )
  end
end

@testset "graded native slot keeps the Hamiltonian branch" begin
  (; Δ, χ, F) = NATIVE_KERR
  d = 3
  model = native_kerr_model(d, Δ, χ, 0.0, F)
  section = native_ladder_section(model, Δ, χ, 0.0)
  graded = native_bf_arbitrary(model, 4; inverse=section)
  ordinary = native_bf_arbitrary(model, 4)

  @test isempty(graded.channels)
  @test all(norm(Sn) <= 1e-7 for Sn in graded.S)
  @test all(isapprox(graded.H[k], ordinary.H[k]; atol=1e-8, rtol=1e-8) for k in 1:5)
end

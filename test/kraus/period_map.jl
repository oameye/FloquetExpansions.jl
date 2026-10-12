using Test
using LinearAlgebra: Hermitian, I, SymTridiagonal, diagm, eigen, eigvals, opnorm
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using FloquetExpansions

const SQA = SecondQuantizedAlgebra

@variables ω::Real t::Real Δ::Real γ::Real

const qubit = NLevelSpace(:period_map_qubit, 2)
σ(i, j) = Transition(qubit, :σ, i, j)
const σz = σ(1, 1) - σ(2, 2)
const σx = σ(1, 2) + σ(2, 1)
const σy = -im * σ(1, 2) + im * σ(2, 1)

include(joinpath(@__DIR__, "..", "helpers", "shared.jl"))

function harmonic_matrices(op, d::Int, substitutions)
  P = harmonics(FloquetExpansions.qadd(op), ω, t)
  return [(m, tomatrix(P[m], d, substitutions)) for m in keys(P)]
end

function no_jump_frame(Φ, algorithm, N)
  lowering = first(FloquetExpansions.period_map_lowering(Φ.generator, N))
  zero_element = zero(FloquetExpansions.algebra_identity(lowering.algebra))
  Q = FloquetExpansions.no_jump_generator(
    FloquetExpansions.lower_harmonics(lowering, Φ.hamiltonian),
    FloquetExpansions.lowered_channels(lowering, Φ),
    zero_element,
  )
  return FloquetExpansions.no_jump_frame(algorithm, Q, N)
end

evaluate(harmonics, phase) = sum(M .* cis(-m * phase) for (m, M) in harmonics)

function lindblad_rate(ρ, Hh, jumps, phase)
  H = evaluate(Hh, phase)
  result = -im * (H * ρ - ρ * H)
  for (rate, Jh) in jumps
    J = evaluate(Jh, phase)
    result += rate * (J * ρ * J' - (J' * J * ρ + ρ * J' * J) / 2)
  end
  return result
end

function exact_period_action(ρ, H, channels, d, value; steps=1200)
  Hh = harmonic_matrices(H, d, Dict())
  jumps = [(Float64(rate), harmonic_matrices(J, d, Dict())) for (J, rate) in channels]
  h = 2π / value / steps
  for k in 0:(steps - 1)
    phase = value * k * h
    k1 = lindblad_rate(ρ, Hh, jumps, phase)
    k2 = lindblad_rate(ρ + h / 2 * k1, Hh, jumps, phase + value * h / 2)
    k3 = lindblad_rate(ρ + h / 2 * k2, Hh, jumps, phase + value * h / 2)
    k4 = lindblad_rate(ρ + h * k3, Hh, jumps, phase + value * h)
    ρ = ρ + h / 6 * (k1 + 2k2 + 2k3 + k4)
  end
  return ρ
end

function map_action(Φ, ρ, d, value)
  substitutions = Dict(ω => value)
  result = zero(ρ)
  for (left, right, c) in terms(superoperator(Φ))
    result +=
      numeric(c, substitutions) *
      (tomatrix(left, d, substitutions) * ρ * tomatrix(right, d, substitutions))
  end
  return result
end

function units(d, block)
  return [(E=zeros(ComplexF64, d, d); E[i, j]=1; E) for i in 1:block, j in 1:block]
end

function period_map_error(Φ, H, channels, value; d=2, block=d)
  return maximum(units(d, block)) do ρ
    difference = map_action(Φ, ρ, d, value) - exact_period_action(ρ, H, channels, d, value)
    return opnorm(difference[1:block, 1:block])
  end
end

convergence_mismatch(error, N) = abs(error(20.0) / error(40.0) / 2^(N + 1) - 1)

sandwich(A, B) = kron(transpose(B), A)

function liouvillian_matrix(L::Liouvillian, substitutions)
  S = zeros(ComplexF64, 4, 4)
  for (left, right, c) in terms(L)
    S .+=
      numeric(c, substitutions) .*
      sandwich(tomatrix(left, 2, substitutions), tomatrix(right, 2, substitutions))
  end
  return S
end

map_matrix(Φ, substitutions) = liouvillian_matrix(superoperator(Φ), substitutions)

function choi_eigenvalues(Φ, d, value)
  choi = zeros(ComplexF64, d^2, d^2)
  for ρ in units(d, d)
    choi .+= kron(ρ, map_action(Φ, ρ, d, value))
  end
  return eigvals(Hermitian(choi))
end

pairs_to_jumps(channels) = Tuple(jump(J, rate) for (J, rate) in channels)

const drive = (1 // 2) * σz + cos(ω * t) * σx
const decay = ((σ(1, 2), 1 // 5), ((1 // 3) * σz, 1 // 1))
const rotating = ((σz + cos(ω * t) * σy + sin(ω * t) * σx, 1 // 4),)

@testset "period map agrees with the exact one-period map through its order" begin
  for (H, channels) in ((drive, ()), (drive, decay), ((1 // 3) * σz, rotating)), N in 1:3
    Φ = floquet_expansion(
      period_map(H, ω, t; channels=pairs_to_jumps(channels)), VanVleck(), N
    )
    @test convergence_mismatch(value -> period_map_error(Φ, H, channels, value), N) < 0.1
  end
end

@testset "period map is completely positive and trace preserving through its order" begin
  Φ = floquet_expansion(
    period_map(drive, ω, t; channels=pairs_to_jumps(decay)), VanVleck(), 3
  )
  for value in (1.0, 2.0, 20.0)
    @test minimum(choi_eigenvalues(Φ, 2, value)) > -1e-10
  end
  defect(value) = opnorm(tomatrix(trace_defect(Φ), 2, Dict(ω => value)))
  @test convergence_mismatch(defect, 3) < 0.1
end

@testset "Hori–Deprit and Bloch/Feshbach give the same period map" begin
  channels = pairs_to_jumps(decay)
  Φ = period_map(drive, ω, t; channels)
  hd = floquet_expansion(Φ, VanVleck(), 3)
  bf = floquet_expansion(Φ, VanVleck(; algorithm=BlochFeshbach()), 3)
  term_table(Φ) = Dict((left, right) => c for (left, right, c) in terms(superoperator(Φ)))
  @test isequal(term_table(hd), term_table(bf))

  H = drive + (1 // 3) * cos(2 * ω * t) * σy
  Φ = period_map(H, ω, t; channels=pairs_to_jumps(rotating))
  hd, bf = (no_jump_frame(Φ, algorithm, 6) for algorithm in (HoriDeprit(), BlochFeshbach()))
  @test all(iszero(a - b) for (a, b) in zip(hd.effective, bf.effective))
  for field in (:kick, :inverse), (a, b) in zip(getfield(hd, field), getfield(bf, field))
    @test all(iszero, values(FloquetExpansions.harmonic_combine(a, b, 1, -1)))
  end
end

function grade_coefficients(f, degree)
  nodes = [cospi((2j + 1) / (2degree + 2)) / 2 for j in 0:degree]
  W = inv([ε^k for ε in nodes, k in 0:degree])
  samples = [f(ε) for ε in nodes]
  return [sum(W[k + 1, j] * samples[j] for j in eachindex(nodes)) for k in 0:degree]
end

truncated_product(A, B, N) = [sum(A[i + 1] * B[k - i + 1] for i in 0:k) for k in 0:N]

function truncated_exp(A, N)
  result = [k == 0 ? Matrix{ComplexF64}(I, size(A[1])) : zero(A[1]) for k in 0:N]
  power = result
  for p in 1:N
    power = truncated_product(power, A, N) ./ p
    result = result .+ power
  end
  return result
end

@testset "period map recovers the effective generator after the micromotion" begin
  for (gauge, H, channels) in
      ((VanVleck(), (1 // 3) * σz, rotating), (GKSLNormalForm(), drive, decay[1:1]))
    N = 3
    Φ = floquet_expansion(
      period_map(H, ω, t; channels=pairs_to_jumps(channels)), VanVleck(), N
    )
    x = floquet_expansion(H, ω, t, gauge, N; channels=pairs_to_jumps(channels))
    at(ε) = Dict(ω => 1 / ε, t => 0)
    degree = 2N + 3
    K = grade_coefficients(ε -> liouvillian_matrix(micromotion(x)(t), at(ε)), degree)
    TL = grade_coefficients(
      ε -> 2π * ε * liouvillian_matrix(effective_generator(x), at(ε)), degree
    )
    S = grade_coefficients(ε -> map_matrix(Φ, at(ε)), degree)
    conjugate = truncated_exp(K[1:(N + 1)], N)
    inverse = truncated_exp(-K[1:(N + 1)], N)
    reference = truncated_product(
      truncated_product(conjugate, truncated_exp(TL[1:(N + 1)], N), N), inverse, N
    )
    for k in 0:N
      @test S[k + 1] ≈ reference[k + 1] atol = 1e-9
    end
  end
end

@testset "period map with symbolic parameters equals the numeric map at substituted points" begin
  symbolic = floquet_expansion(
    period_map(Δ * σz + cos(ω * t) * σx, ω, t; channels=(jump(σ(1, 2), γ),)), VanVleck(), 2
  )
  exact = floquet_expansion(
    period_map(drive, ω, t; channels=(jump(σ(1, 2), 1 // 5),)), VanVleck(), 2
  )
  point = Dict(ω => 9.0, Δ => 1 / 2, γ => 1 / 5)
  @test map_matrix(symbolic, point) ≈ map_matrix(exact, point) atol = 1e-12
end

@testset "period map of driven Kerr with loss matches the truncated dynamics" begin
  fock = FockSpace(:period_map_cavity)
  a = Destroy(fock, :a)
  H = (1 // 2) * a' * a + (3 // 10) * a' * a * a' * a + (2 // 5) * (a + a') * cos(ω * t)
  channels = ((a, 4 // 5),)
  N = 2
  Φ = floquet_expansion(
    period_map(H, ω, t; channels=pairs_to_jumps(channels)), VanVleck(), N
  )
  error(value) = period_map_error(Φ, H, channels, value; d=14, block=3)
  @test convergence_mismatch(error, N) < 0.15
end

function gauss_legendre(n)
  decomposition = eigen(SymTridiagonal(zeros(n), [k / sqrt(4k^2 - 1) for k in 1:(n - 1)]))
  return (decomposition.values .+ 1) ./ 2, decomposition.vectors[1, :] .^ 2
end

function simplex_rule(q, n)
  nodes, weights = gauss_legendre(n)
  rule = [(Float64[], 1.0)]
  for _ in 1:q
    rule = [
      begin
        start = isempty(u) ? 0.0 : u[end]
        ([u; start + (1 - start) * x], w * (1 - start) * v)
      end for (u, w) in rule for (x, v) in zip(nodes, weights)
    ]
  end
  return rule
end

function amplitude_action(Φ, ρ, value)
  result = zero(ρ)
  for K in jump_amplitudes(Φ)
    w = numeric(weight(K), Dict(ω => value))
    for (u, v) in simplex_rule(length(jump_times(K)), 24)
      M = tomatrix(amplitude(K), 2, Dict(ω => value, (jump_times(K) .=> u)...))
      result += (w * v) * (M * ρ * M')
    end
  end
  return result
end

@testset "jump amplitudes integrate to the period map" begin
  Φ = floquet_expansion(
    period_map(drive, ω, t; channels=pairs_to_jumps(rotating)), VanVleck(), 2
  )
  for ρ in units(2, 2)
    @test amplitude_action(Φ, ρ, 3.0) ≈ map_action(Φ, ρ, 2, 3.0) atol = 1e-12
  end
end

const mixing = ((3 // 5, 4 // 5), (-4 // 5, 3 // 5))

@testset "period map is covariant under its representation freedoms" begin
  J = (σ(1, 2), (1 // 2) * σz)
  mixed = Tuple(jump(a * J[1] + b * J[2], 1 // 3) for (a, b) in mixing)
  plain = Tuple(jump(L, 1 // 3) for L in J)
  for value in (2.0, 9.0)
    point = Dict(ω => value)
    @test map_matrix(
      floquet_expansion(period_map(drive, ω, t; channels=mixed), VanVleck(), 2), point
    ) ≈ map_matrix(
      floquet_expansion(period_map(drive, ω, t; channels=plain), VanVleck(), 2), point
    ) atol = 1e-12
  end
  c = 1 // 2
  shifted_H = drive + (im * c / 2) * (σ(2, 1) - σ(1, 2)) * (1 // 3)
  N = 2
  shifted = floquet_expansion(
    period_map(shifted_H, ω, t; channels=(jump(σ(1, 2) + c, 1 // 3),)), VanVleck(), N
  )
  original = floquet_expansion(
    period_map(drive, ω, t; channels=(jump(σ(1, 2), 1 // 3),)), VanVleck(), N
  )
  difference(value) =
    opnorm(map_matrix(shifted, Dict(ω => value)) - map_matrix(original, Dict(ω => value)))
  @test difference(20.0) / difference(40.0) > 0.9 * 2^(N + 1)
end

@testset "period map rejects what it cannot expand" begin
  @test_throws ArgumentError floquet_expansion(period_map(drive, ω, t), VanVleck(), 0)
  @test_throws ArgumentError floquet_expansion(period_map(drive, ω, t), GKSLNormalForm(), 1)
  @test_throws ArgumentError period_map(
    drive, ω, t; channels=(jump(σ(1, 2), (1 // 5) * (1 + cos(ω * t))),)
  )
  @test_throws ArgumentError floquet_expansion(period_map(ω * σz, ω, t), VanVleck(), 1)
  @test_throws ArgumentError floquet_expansion(
    period_map(drive, ω, t; channels=(jump(σ(1, 2), 0.2),)), VanVleck(), 1
  )
end

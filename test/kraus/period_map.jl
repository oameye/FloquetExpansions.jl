using Test
using LinearAlgebra: Hermitian, I, diagm, eigvals, opnorm
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

function numeric(c::SQA.CNum, substitutions)
  value(x) =
    ComplexF64(Symbolics.value(Symbolics.substitute(x, substitutions; fold=Val(true))))
  z = SQA.to_num(c)
  return value(real(z)) + im * value(imag(z))
end

function operator_matrix(op::SQA.Op, d::Int)
  lowering = diagm(1 => sqrt.(1.0:(d - 1)))
  SQA.is_destroy(op) && return complex(lowering)
  SQA.is_create(op) && return complex(lowering')
  E = zeros(ComplexF64, d, d)
  E[op.l1, op.l2] = 1
  return E
end

function dense(q, d::Int, substitutions)
  M = zeros(ComplexF64, d, d)
  for (term, c) in FloquetExpansions.qadd(q)
    product = Matrix{ComplexF64}(I, d, d)
    for op in term.ops
      product = product * operator_matrix(op, d)
    end
    M .+= numeric(c, substitutions) .* product
  end
  return M
end

function harmonic_matrices(op, d::Int, substitutions)
  P = harmonics(FloquetExpansions.qadd(op), ω, t)
  return [(m, dense(P[m], d, substitutions)) for m in keys(P)]
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
      (dense(left, d, substitutions) * ρ * dense(right, d, substitutions))
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
      sandwich(dense(left, 2, substitutions), dense(right, 2, substitutions))
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
    Φ = period_map(H, ω, t, N; channels=pairs_to_jumps(channels))
    @test convergence_mismatch(value -> period_map_error(Φ, H, channels, value), N) < 0.1
  end
end

@testset "period map is completely positive and trace preserving through its order" begin
  Φ = period_map(drive, ω, t, 3; channels=pairs_to_jumps(decay))
  for value in (1.0, 2.0, 20.0)
    @test minimum(choi_eigenvalues(Φ, 2, value)) > -1e-10
  end
  defect(value) = opnorm(dense(trace_defect(Φ), 2, Dict(ω => value)))
  @test convergence_mismatch(defect, 3) < 0.1
end

@testset "Hori–Deprit and Bloch/Feshbach give the same period map" begin
  channels = pairs_to_jumps(decay)
  hd = period_map(drive, ω, t, 3; channels)
  bf = period_map(drive, ω, t, 3; channels, algorithm=BlochFeshbach())
  term_table(Φ) = Dict((left, right) => c for (left, right, c) in terms(superoperator(Φ)))
  @test isequal(term_table(hd), term_table(bf))
end

@testset "period map recovers the effective generator after the micromotion" begin
  for (gauge, H, channels) in
      ((VanVleck(), (1 // 3) * σz, rotating), (GKSLNormalForm(), drive, decay[1:1]))
    N = 3
    Φ = period_map(H, ω, t, N; channels=pairs_to_jumps(channels))
    x = floquet_expansion(H, ω, t, gauge, N; channels=pairs_to_jumps(channels))
    function error(value)
      substitutions = Dict(ω => value, t => 0)
      K = liouvillian_matrix(micromotion(x)(t), substitutions)
      L = liouvillian_matrix(effective_generator(x), substitutions)
      reference = exp(K) * exp(2π / value * L) * exp(-K)
      return opnorm(map_matrix(Φ, substitutions) - reference)
    end
    @test convergence_mismatch(error, N) < 0.1
  end
end

@testset "period map with symbolic parameters equals the numeric map at substituted points" begin
  symbolic = period_map(Δ * σz + cos(ω * t) * σx, ω, t, 2; channels=(jump(σ(1, 2), γ),))
  exact = period_map(drive, ω, t, 2; channels=(jump(σ(1, 2), 1 // 5),))
  point = Dict(ω => 9.0, Δ => 1 / 2, γ => 1 / 5)
  @test map_matrix(symbolic, point) ≈ map_matrix(exact, point) atol = 1e-12
end

@testset "period map of driven Kerr with loss matches the truncated dynamics" begin
  fock = FockSpace(:period_map_cavity)
  a = Destroy(fock, :a)
  H = (1 // 2) * a' * a + (3 // 10) * a' * a * a' * a + (2 // 5) * (a + a') * cos(ω * t)
  channels = ((a, 4 // 5),)
  N = 2
  Φ = period_map(H, ω, t, N; channels=pairs_to_jumps(channels))
  error(value) = period_map_error(Φ, H, channels, value; d=14, block=3)
  @test convergence_mismatch(error, N) < 0.15
end

@testset "period map rejects what it cannot expand" begin
  @test_throws ArgumentError period_map(drive, ω, t, 0)
  @test_throws ArgumentError period_map(
    drive, ω, t, 1; channels=(jump(σ(1, 2), (1 // 5) * (1 + cos(ω * t))),)
  )
  @test_throws ArgumentError period_map(ω * σz, ω, t, 1)
  @test_throws ArgumentError period_map(drive, ω, t, 1; channels=(jump(σ(1, 2), 0.2),))
end

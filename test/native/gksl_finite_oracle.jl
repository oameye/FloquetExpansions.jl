using Test
using LinearAlgebra: Hermitian, I, eigmin, eigvals, norm
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using FloquetExpansions

const SQA = SecondQuantizedAlgebra

@variables ω::Real t::Real

include(joinpath(@__DIR__, "..", "helpers", "shared.jl"))

function dense_rate(rate, value::Real)
  return real(ComplexF64(SQA.to_complex(SQA.substitute(rate, Dict(ω => value)))))
end

function hamiltonian_superoperator(H::AbstractMatrix)
  Id = Matrix{ComplexF64}(I, size(H)...)
  return -im * (kron(Id, H) - kron(transpose(H), Id))
end

function dissipator_superoperator(L::AbstractMatrix)
  Id = Matrix{ComplexF64}(I, size(L)...)
  LL = L' * L
  return kron(conj(L), L) - kron(Id, LL) / 2 - kron(transpose(LL), Id) / 2
end

function effective_superoperator(expansion, d::Int, value::Real)
  S = hamiltonian_superoperator(tomatrix(hamiltonian(expansion), d, Dict(ω => value)))
  for (channel, rate) in channels(expansion)
    S +=
      dense_rate(rate, value) *
      dissipator_superoperator(tomatrix(SQA.jump(channel), d, Dict(ω => value)))
  end
  return S
end

function numeric_kossakowski(expansion, value::Real)
  return [
    ComplexF64(SQA.to_complex(SQA.substitute(c, Dict(ω => value)))) for
    c in kossakowski(expansion)
  ]
end

function rk4_monodromy(rhs, period::Float64, steps::Int, dim::Int)
  U = Matrix{ComplexF64}(I, dim, dim)
  h = period / steps
  for k in 0:(steps - 1)
    s = k * h
    k1 = rhs(s) * U
    k2 = rhs(s + h / 2) * (U + h / 2 * k1)
    k3 = rhs(s + h / 2) * (U + h / 2 * k2)
    k4 = rhs(s + h) * (U + h * k3)
    U += h / 6 * (k1 + 2k2 + 2k3 + k4)
  end
  return U
end

function floquet_exponents(static, drive, value::Real)
  period = 2π / value
  rhs(s) = static + cos(value * s) * drive
  return log.(eigvals(rk4_monodromy(rhs, period, 4000, size(static, 1)))) ./ period
end

function spectral_distance(exact, approximate)
  remaining = collect(approximate)
  worst = 0.0
  for λ in exact
    index = argmin(abs.(remaining .- λ))
    worst = max(worst, abs(remaining[index] - λ))
    deleteat!(remaining, index)
  end
  return worst
end

function loglog_slope(values, errors)
  return -(log(errors[end]) - log(errors[1])) / (log(values[end]) - log(values[1]))
end

h = NLevelSpace(:gksl_oracle, 3)
σ(i, j) = Transition(h, :σ, i, j)
H0 = (1 // 2) * σ(2, 2) + (6 // 5) * σ(3, 3)
X = σ(1, 2) + σ(2, 1) + (3 // 4) * (σ(2, 3) + σ(3, 2))
H = H0 + (2 // 5) * X * cos(ω * t)
decay = ((σ(1, 2), 3 // 10), (σ(2, 3), 1 // 5), (σ(2, 2), 1 // 10))

function lindbladian(H, channels)
  return foldl(channels; init=hamiltonian_action(H)) do L, (J, rate)
    return L + rate * dissipator(J)
  end
end

function lab_static(d)
  S = hamiltonian_superoperator(tomatrix(H0, d, Dict(ω => 1.0)))
  for (J, rate) in decay
    S += Float64(rate) * dissipator_superoperator(tomatrix(J, d, Dict(ω => 1.0)))
  end
  return S
end
lab_drive(d) = hamiltonian_superoperator(tomatrix((2 // 5) * X, d, Dict(ω => 1.0)))

const FREQUENCIES = [20.0, 40.0, 80.0]

@testset "GKSLNormalForm spectrum matches the exact Floquet exponents" begin
  static, drive = lab_static(3), lab_drive(3)
  exact = [floquet_exponents(static, drive, value) for value in FREQUENCIES]
  for order in (2, 3, 4)
    expansion = floquet_expansion(lindbladian(H, decay), ω, t, GKSLNormalForm(), order)
    errors = [
      spectral_distance(exact[k], eigvals(effective_superoperator(expansion, 3, value))) for
      (k, value) in enumerate(FREQUENCIES)
    ]
    @test issorted(errors; rev=true)
    @test loglog_slope(FREQUENCIES, errors) >= order - 0.1
  end
end

@testset "GKSLNormalForm driven Kerr matches the truncated Floquet exponents" begin
  mode = FockSpace(:gksl_oracle_kerr)
  a = Destroy(mode, :a)
  K0 = (1 // 2) * a' * a + (3 // 10) * a' * a * a' * a
  V = (2 // 5) * (a + a')
  L = hamiltonian_action(K0 + V * cos(ω * t)) + (4 // 5) * dissipator(a)
  expansions = [floquet_expansion(L, ω, t, GKSLNormalForm(), order) for order in (4, 5)]
  for d in (6, 8)
    static =
      hamiltonian_superoperator(tomatrix(K0, d)) +
      (4 / 5) * dissipator_superoperator(tomatrix(a, d))
    drive = hamiltonian_superoperator(tomatrix(V, d))
    slowest = [
      sort(floquet_exponents(static, drive, value); by=λ -> -real(λ))[1:6] for
      value in FREQUENCIES
    ]
    for (order, expansion) in zip((4, 5), expansions)
      errors = [
        spectral_distance(
          slowest[k], eigvals(effective_superoperator(expansion, d, value))
        ) for (k, value) in enumerate(FREQUENCIES)
      ]
      @test issorted(errors; rev=true)
      @test loglog_slope(FREQUENCIES, errors) >= order - 0.1
    end
  end
end

@testset "GKSLNormalForm Kossakowski matrix is positive semidefinite at finite frequency" begin
  qubit = NLevelSpace(:gksl_oracle_qubit, 2)
  τ(i, j) = Transition(qubit, :τ, i, j)
  Hq = (1 // 2) * τ(2, 2) + (2 // 5) * (τ(1, 2) + τ(2, 1)) * cos(ω * t)
  models = (
    (lindbladian(H, decay), (3, 4)),
    (lindbladian(Hq, ((τ(1, 2), 3 // 10), (τ(2, 2), 1 // 10))), (3, 4, 5, 6)),
  )
  for (generator, orders) in models, order in orders
    expansion = floquet_expansion(generator, ω, t, GKSLNormalForm(), order)
    for value in (10.0, 20.0, 40.0, 80.0)
      C = numeric_kossakowski(expansion, value)
      @test norm(C - C') <= 1.0e-12 * max(1, norm(C))
      @test eigmin(Hermitian(C)) >= -1.0e-12
    end
  end
end

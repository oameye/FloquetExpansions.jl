using Test
using LinearAlgebra: Hermitian, I, eigmin, eigvals, norm
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using FloquetExpansions

const SQA = SecondQuantizedAlgebra

@variables ω::Real t::Real

function dense_operator(q::SQA.QAdd, d::Int, value::Real)
  M = zeros(ComplexF64, d, d)
  for (term, coeff) in SQA.substitute(q, Dict(ω => value))
    z = SQA.to_num(coeff)
    c = complex(Float64(Symbolics.value(real(z))), Float64(Symbolics.value(imag(z))))
    T = Matrix{ComplexF64}(I, d, d)
    for o in term.ops
      E = zeros(ComplexF64, d, d)
      E[o.l1, o.l2] = 1
      T = T * E
    end
    M .+= c .* T
  end
  return M
end

function dense_operator(op::SQA.Op, d::Int, ::Real)
  E = zeros(ComplexF64, d, d)
  E[op.l1, op.l2] = 1
  return E
end

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
  S = hamiltonian_superoperator(dense_operator(hamiltonian(expansion), d, value))
  for c in channels(expansion)
    S +=
      dense_rate(c.rate, value) *
      dissipator_superoperator(dense_operator(c.operator, d, value))
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
decay = (jump(σ(1, 2), 3 // 10), jump(σ(2, 3), 1 // 5), jump(σ(2, 2), 1 // 10))

function lab_static(d)
  S = hamiltonian_superoperator(dense_operator(H0, d, 1.0))
  for c in decay
    S +=
      dense_rate(c.rate, 1.0) * dissipator_superoperator(dense_operator(c.operator, d, 1.0))
  end
  return S
end
lab_drive(d) = hamiltonian_superoperator(dense_operator((2 // 5) * X, d, 1.0))

const FREQUENCIES = [20.0, 40.0, 80.0]

@testset "GKSLNormalForm spectrum matches the exact Floquet exponents" begin
  static, drive = lab_static(3), lab_drive(3)
  exact = [floquet_exponents(static, drive, value) for value in FREQUENCIES]
  for order in (2, 3, 4)
    expansion = floquet_expansion(H, ω, t, GKSLNormalForm(), order; channels=decay)
    errors = [
      spectral_distance(exact[k], eigvals(effective_superoperator(expansion, 3, value))) for
      (k, value) in enumerate(FREQUENCIES)
    ]
    @test issorted(errors; rev=true)
    @test loglog_slope(FREQUENCIES, errors) >= order - 0.1
  end
end

@testset "GKSLNormalForm Kossakowski matrix is positive semidefinite at finite frequency" begin
  qubit = NLevelSpace(:gksl_oracle_qubit, 2)
  τ(i, j) = Transition(qubit, :τ, i, j)
  Hq = (1 // 2) * τ(2, 2) + (2 // 5) * (τ(1, 2) + τ(2, 1)) * cos(ω * t)
  models = (
    (H, decay, (3, 4)), (Hq, (jump(τ(1, 2), 3 // 10), jump(τ(2, 2), 1 // 10)), (3, 4, 5, 6))
  )
  for (generator, model_channels, orders) in models, order in orders
    expansion = floquet_expansion(
      generator, ω, t, GKSLNormalForm(), order; channels=model_channels
    )
    for value in (10.0, 20.0, 40.0, 80.0)
      C = numeric_kossakowski(expansion, value)
      @test norm(C - C') <= 1.0e-12 * max(1, norm(C))
      @test eigmin(Hermitian(C)) >= -1.0e-12
    end
  end
end

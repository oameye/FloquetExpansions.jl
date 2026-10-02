using Test
using LinearAlgebra: norm
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions
const Q = Complex{Rational{BigInt}}
const R = Rational{BigInt}

# Driven Kerr oscillator with loss in the unnormalized Fock basis |n) = a†^n |0>, Hilbert metric n!.
function exact_kerr(d)
  g = R[factorial(big(n)) for n in 0:(d - 1)]
  rep = FE.ExactLiouvilleRepresentation{Q}(g)
  a = zeros(Q, d, d)
  for n in 1:(d - 1)
    a[n, n + 1] = n
  end
  ad = FE.hilbert_adjoint(g, a)
  n = ad * a
  c = (rep.inverse_frame * vec(a))[2:end]
  Z = zeros(Q, d^2 - 1, d^2 - 1)
  F = 2 // 5 + im // 5
  L = Dict{Int,Matrix{Q}}(
    0 => FE.native_gksl(rep, n / 2 + (3 // 10) * n * n, (4 // 5) * c * c'),
    1 => FE.native_gksl(rep, F * (a + ad), Z),
    -1 => FE.native_gksl(rep, conj(F) * (a + ad), Z),
  )
  return rep, L, reshape(c, :, 1)
end

@testset "charge grading detects the Fock number and solves the charged sectors exactly" begin
  d = 4
  rep, L, leading = exact_kerr(d)
  inverse = FE.ChargeGradedInverse(rep, L[0], leading, 1e-9)
  @test inverse.weights == [0 1 2 3]
  @test inverse.singular == [[0]]
  @test !isempty(inverse.regular)

  # A charged dark target: the regular gauge cancels it exactly and leaves the neutral part.
  dim = d^2 - 1
  X = zeros(Q, d, d)
  X[1, 3] = 1
  Y = zeros(Q, d, d)
  Y[2, 3] = 1 // 2
  cx = (rep.inverse_frame * vec(X))[2:end]
  cy = (rep.inverse_frame * vec(Y))[2:end]
  C = cx * cy' + cy * cx' + cy * cy'
  S = FE.regular_gauge(inverse, rep, L[0], C, zeros(Q, dim, dim), leading, 1e-9)
  target = FE.native_gksl(
    rep, zeros(Q, d, d), FE.native_dark_target(rep, C, zeros(Q, dim, dim), leading, 1e-9)
  )
  image = L[0] * S - S * L[0]
  charged = [
    (o, i) for o in 1:(d ^ 2), i in 1:(d ^ 2) if
    !iszero(FE.superoperator_label(inverse.weights, d, o, i))
  ]
  @test all(image[o, i] == -target[o, i] for (o, i) in charged)
  @test all(iszero(S[o, i]) for o in 1:(d ^ 2), i in 1:(d ^ 2) if (o, i) ∉ charged)
end

@testset "charge grading falls back to one sector without a symmetry" begin
  rep = FE.ExactLiouvilleRepresentation{Q}(2)
  L0 = FE.native_gksl(rep, Q[0 1; 1 0], zeros(Q, 3, 3))
  inverse = FE.ChargeGradedInverse(rep, L0, zeros(Q, 3, 0), 1e-9)
  @test size(inverse.weights, 1) == 0
  @test isempty(inverse.regular)
  @test length(inverse.directions) == length(rep.gauge)
end

@testset "charge-graded recurrence on exact driven Kerr through order one" begin
  rep, L, leading = exact_kerr(4)
  inverse = FE.ChargeGradedInverse(rep, L[0], leading, 1e-9)
  bf = FE.native_recurrence(
    FE.BlochFeshbach(), rep, inverse, L, leading, Q[4 // 5], 1, 1e-8
  )
  hd = FE.native_recurrence(FE.HoriDeprit(), rep, inverse, L, leading, Q[4 // 5], 1, 1e-8)
  @test bf.E == hd.E
end

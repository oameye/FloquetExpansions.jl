using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra
using Symbolics: Symbolics
using LinearAlgebra: exp, opnorm
using Random: MersenneTwister, randn

include(joinpath(@__DIR__, "helpers", "shared.jl"))

const D = 3
const HN = NLevelSpace(:atom, D)

function tomatrices(X::PeriodicGenerator, d::Int)
  return Dict{Int,Matrix{ComplexF64}}(l => tomatrix(X[l], d) for l in keys(X))
end

function frechet_exp(A::Matrix{ComplexF64}, E::Matrix{ComplexF64})
  d = size(A, 1)
  B = zeros(ComplexF64, 2d, 2d)
  B[1:d, 1:d] = A
  B[1:d, (d + 1):(2d)] = E
  B[(d + 1):(2d), (d + 1):(2d)] = A
  X = exp(B)
  return X[1:d, 1:d], X[1:d, (d + 1):(2d)]
end

function residual(vv, H, wd, d; nt=24)
  Heff = tomatrix(effective_generator(vv), d)
  Ks = tomatrices(micromotion(vv), d)
  Hs = tomatrices(H, d)
  worst = 0.0
  for k in 0:(nt - 1)
    t = 2π * k / (nt * wd)
    Kt = zeros(ComplexF64, d, d)
    Kdt = zeros(ComplexF64, d, d)
    Ht = zeros(ComplexF64, d, d)
    for (l, Kl) in Ks
      ph = cis(-l * wd * t)
      Kt .+= Kl .* ph
      Kdt .+= Kl .* (-im * l * wd * ph)
    end
    for (m, Hm) in Hs
      Ht .+= Hm .* cis(-m * wd * t)
    end
    A = -im .* Kt
    expA, L = frechet_exp(A, -im .* Kdt)
    expmA = exp(-A)
    worst = max(worst, opnorm(expmA * Ht * expA .- im .* (expmA * L) .- Heff))
  end
  return worst
end

function fitslope(xs, ys)
  lx, ly = log.(xs), log.(ys)
  mx, my = sum(lx) / length(lx), sum(ly) / length(ly)
  return sum((lx .- mx) .* (ly .- my)) / sum((lx .- mx) .^ 2)
end

@testset "the truncated factorization solves eq:defining0 to O(wd^-(N-1))" begin
  wds = [20.0, 40.0, 80.0, 160.0]

  for N in 2:5
    errs = [
      let H = random_drive(MersenneTwister(0xF10), HN, D, 2, wd)
        residual(floquet_expansion(H, VanVleck(), N), H, wd, D)
      end for wd in wds
    ]
    @test fitslope(wds, errs) ≈ -(N - 1) atol = 0.1
    @test issorted(errs; rev=true)
  end

  rwa = [
    let H = random_drive(MersenneTwister(0xF10), HN, D, 2, wd)
      residual(floquet_expansion(H, VanVleck(), 1), H, wd, D)
    end for wd in wds
  ]
  @test fitslope(wds, rwa) ≈ 0 atol = 0.1
end

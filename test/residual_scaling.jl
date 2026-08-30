using Test
using FloquetExpansions
import SecondQuantizedAlgebra as SQA
using Symbolics: Symbolics
using LinearAlgebra: I, exp, opnorm
using Random: MersenneTwister, randn

# THE primary gate (DESIGN §8). Everything else either stops at order 2 (the spec's closed forms)
# or checks the triangle against its own definition. This checks K and Heff JOINTLY, at every
# order, against eq:defining0 itself:
#
#     Heff = e^{iK} H_S e^{-iK} - i e^{iK} d_t e^{-iK}
#
# Truncating at order N must leave a residual of O(wd^-(N-1)) — one order less than N because
# eq:defining0 contains d_t, which costs an order, and the omitted Kdot^(N) is an order-(N-1)
# object. Needs no oracle, no license, no transcription.
#
# Verified SHARP against two corruptions of the engine, each of which flattens every slope to
# -1.00: the Deprit weight i/(j+1) -> i/j on the Kdot family, and the k-range n-j+1 -> n-j.

const D = 3
const HN = NLevelSpace(:atom, D)

# `Op`'s docstring fixes the field layout: a transition packs (l1, l2, g, nlev) as
# (i, j, ground state, number of levels), so sigma_ij is the matrix unit E_ij.
function tomatrix(q::SQA.QAdd, d::Int)
  M = zeros(ComplexF64, d, d)
  for (term, coeff) in q
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

function tomatrices(X::PeriodicOperator, d::Int)
  return Dict{Int,Matrix{ComplexF64}}(l => tomatrix(X[l], d) for l in keys(X))
end

# exp'(A)[E] from expm([A E; 0 A]) = [e^A  L; 0  e^A], so d/dt exp(-iK(t)) needs no extra package.
function frechet_exp(A::Matrix{ComplexF64}, E::Matrix{ComplexF64})
  d = size(A, 1)
  B = zeros(ComplexF64, 2d, 2d)
  B[1:d, 1:d] = A
  B[1:d, (d + 1):(2d)] = E
  B[(d + 1):(2d), (d + 1):(2d)] = A
  X = exp(B)
  return X[1:d, 1:d], X[1:d, (d + 1):(2d)]
end

function random_drive(rng, d, M, wd)
  randop() =
    sum(complex(randn(rng), randn(rng)) * Transition(HN, :σ, i, j) for i in 1:d, j in 1:d)
  comps = Dict{Int,SQA.QAdd}()
  h0 = randop()
  comps[0] = h0 + adjoint(h0)
  for m in 1:M
    hm = randop()
    comps[m] = hm
    comps[-m] = adjoint(hm)                    # eq:fourierH
  end
  return PeriodicOperator(comps, wd)
end

function residual(vv, H, wd, d; nt=24)
  Heff = tomatrix(effective_hamiltonian(vv), d)
  Ks = tomatrices(kick_operator(vv), d)
  Hs = tomatrices(H, d)
  worst = 0.0
  for k in 0:(nt - 1)
    t = 2π * k / (nt * wd)                     # sample one drive period
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
      let H = random_drive(MersenneTwister(0xF10), D, 2, wd)
        residual(floquet_expansion(H, VanVleck(), N), H, wd, D)
      end for wd in wds
    ]
    @test fitslope(wds, errs) ≈ -(N - 1) atol = 0.1
    @test issorted(errs; rev=true)           # monotone, not just a good fit
  end

  # N = 1 is the RWA: dropping every 1/wd correction leaves a residual that does not shrink.
  rwa = [
    let H = random_drive(MersenneTwister(0xF10), D, 2, wd)
      residual(floquet_expansion(H, VanVleck(), 1), H, wd, D)
    end for wd in wds
  ]
  @test fitslope(wds, rwa) ≈ 0 atol = 0.1
end

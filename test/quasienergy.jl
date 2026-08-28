using Test
using FloquetExpansions
import SecondQuantizedAlgebra as SQA
using Symbolics: Symbolics
using LinearAlgebra: I, eigvals
using Random: MersenneTwister, randn

# `QuasienergyOperator` is an inspection view: the expansion never builds a Sambe matrix, so
# nothing else in the suite exercises it and it needs its own test (DESIGN §3).
#
# The real check is against exact Floquet theory. Q's eigenvalues ARE the quasienergies, and
# effective_hamiltonian is precisely what approximates them, so diagonalizing Q and comparing
# spectra validates the sign convention and the whole expansion at once, from a direction no
# other test approaches.

const D = 2
const HN = NLevelSpace(:atom, D)

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

function random_drive(rng, d, M)
  randop() =
    sum(complex(randn(rng), randn(rng)) * Transition(HN, :σ, i, j) for i in 1:d, j in 1:d)
  comps = Dict{Int,SQA.QAdd}()
  h0 = randop()
  comps[0] = h0 + adjoint(h0)
  for m in 1:M
    hm = randop()
    comps[m] = hm
    comps[-m] = adjoint(hm)
  end
  return PeriodicOperator(comps)
end

# assemble the truncated Sambe matrix, blocks laid out in harmonic order
function sambe(Q::QuasienergyOperator, d::Int)
  r = harmonic_range(Q)
  n = length(r)
  M = zeros(ComplexF64, n * d, n * d)
  for (i, m) in enumerate(r), (j, k) in enumerate(r)
    M[((i - 1) * d + 1):(i * d), ((j - 1) * d + 1):(j * d)] = tomatrix(Q[m, k], d)
  end
  return M
end

fold(x, wd) = mod(real(x) + wd / 2, wd) - wd / 2

@testset "indexing is by harmonic, and the diagonal carries -m*wd" begin
  @variables w::Real
  a = Transition(HN, :σ, 1, 2)
  H = PeriodicOperator(0 => 1 * (a' * a), 1 => 1 * a, -1 => 1 * a')
  Q = QuasienergyOperator(H, w, 2)

  @test harmonic_range(Q) == -2:2
  @test size(Q) == (5, 5)
  @test iszero(SQA.simplify(Q[1, 0] - H[1]))           # off-diagonal is H_{m-n}
  @test iszero(SQA.simplify(Q[0, 1] - H[-1]))
  @test iszero(SQA.simplify(Q[2, -2] - H[4]))          # zero, |m-n| exceeds the drive
  @test iszero(SQA.simplify(Q[0, 0] - H[0]))           # m = 0 carries no shift

  # MINUS on the diagonal. The opposite Fourier convention would put +m*wd here, and would
  # also transpose the off-diagonal; getting this wrong is invisible without a spectral test.
  @test iszero(SQA.simplify(Q[1, 1] - (H[0] - w * one(SQA.QAdd))))
  @test iszero(SQA.simplify(Q[-2, -2] - (H[0] + 2w * one(SQA.QAdd))))

  @test_throws BoundsError Q[3, 0]
  @test_throws ArgumentError QuasienergyOperator(H, w, -1)
end

@testset "quasienergies match the effective Hamiltonian's spectrum" begin
  H = random_drive(MersenneTwister(0x51E), D, 1)
  wd = 60.0

  quasienergies(nmax) = sort!(
    unique!(
      round.(
        [fold(z, wd) for z in eigvals(sambe(QuasienergyOperator(H, wd, nmax), D))];
        digits=8,
      ),
    ),
  )

  # The Sambe truncation is not what limits accuracy here: a drive with finitely many harmonics
  # converges in nmax long before the expansion error matters. Measured identical to 8 digits
  # from nmax = 10 to 28, so nmax is not a tuning knob hiding a failure.
  @test quasienergies(10) == quasienergies(20)

  Qs = quasienergies(14)
  errs = Float64[]
  for N in 1:5
    vv = floquet_expansion(H, wd, VanVleck(), N)
    heff = sort!([fold(z, wd) for z in eigvals(tomatrix(effective_hamiltonian(vv), D))])
    push!(errs, maximum(minimum(abs(e - q) for q in Qs) for e in heff))
  end

  # Measured 0.0293, 1.36e-3, 6.73e-6, 1.46e-6, 4.46e-9. Monotone, but the per-order ratio is
  # irregular (21.6, 202, 4.6, 328) because the error is a max over eigenvalues and individual
  # levels can near-cancel at a given order. So assert the envelope, not a per-step rate.
  @test issorted(errs; rev=true)
  @test errs[end] < 1.0e-7
  @test errs[1] / errs[end] > 1.0e5

  # This is what pins the MINUS on the diagonal: rebuilt with +m*wd, the error saturates at
  # 0.056 and stops converging at every order past the first.
end

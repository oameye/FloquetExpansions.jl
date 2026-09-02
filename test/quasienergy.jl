using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra
using Symbolics: Symbolics
using LinearAlgebra: eigvals
using Random: MersenneTwister, randn

include(joinpath(@__DIR__, "helpers", "shared.jl"))

const D = 2
const HN = NLevelSpace(:atom, D)
@variables wd_symbol::Real

function sambe(Q::QuasienergyOperator, d::Int, substitutions::AbstractDict)
  r = harmonic_range(Q)
  n = length(r)
  M = zeros(ComplexF64, n * d, n * d)
  for (i, m) in enumerate(r), (j, k) in enumerate(r)
    M[((i - 1) * d + 1):(i * d), ((j - 1) * d + 1):(j * d)] = tomatrix(
      Q[m, k], d, substitutions
    )
  end
  return M
end

fold(x, wd) = mod(real(x) + wd / 2, wd) - wd / 2

@testset "indexing is by harmonic, and the diagonal carries -m*wd" begin
  @variables w::Real
  a = Transition(HN, :σ, 1, 2)
  H = PeriodicGenerator(Dict(0 => 1 * (a' * a), 1 => 1 * a, -1 => 1 * a'), w)
  Q = QuasienergyOperator(H, 2)

  @test harmonic_range(Q) == -2:2
  @test size(Q) == (5, 5)
  @test iszero(SQA.simplify(Q[1, 0] - H[1]))           # off-diagonal is H_{m-n}
  @test iszero(SQA.simplify(Q[0, 1] - H[-1]))
  @test iszero(SQA.simplify(Q[2, -2] - H[4]))          # zero, |m-n| exceeds the drive
  @test iszero(SQA.simplify(Q[0, 0] - H[0]))

  @test iszero(SQA.simplify(Q[1, 1] - (H[0] - w * one(SQA.QAdd))))
  @test iszero(SQA.simplify(Q[-2, -2] - (H[0] + 2w * one(SQA.QAdd))))

  @test_throws BoundsError Q[3, 0]
  @test_throws ArgumentError QuasienergyOperator(H, -1)
end

@testset "dissipative quasienergy blocks use the Liouvillian convention" begin
  space = NLevelSpace(:dissipative_quasienergy, 2)
  a = Transition(space, :σ, 1, 2)
  static = Liouvillian(zero(SQA.QAdd); channels=(jump(a, 1),))
  driven = hamiltonian_action(a)
  L = PeriodicGenerator(Dict(0 => static, 1 => driven), wd_symbol)
  Q = QuasienergyOperator(L, 1)
  rho = Transition(space, :rho, 1, 1)

  @test Q isa QuasienergyOperator{Liouvillian}
  @test Q[1, 0] == im * driven
  @test iszero(SQA.simplify(Q[1, 1](rho) - (im * static(rho) - wd_symbol * rho)))
end

@testset "quasienergies match the effective Hamiltonian's spectrum" begin
  wd = 60.0
  H = random_drive(MersenneTwister(0x51E), HN, D, 1, wd_symbol)
  substitutions = Dict(wd_symbol => wd)

  quasienergies(nmax) = sort!(
    unique!(
      round.(
        [
          fold(z, wd) for
          z in eigvals(sambe(QuasienergyOperator(H, nmax), D, substitutions))
        ];
        digits=8,
      ),
    ),
  )

  @test quasienergies(10) == quasienergies(20)

  Qs = quasienergies(14)
  errs = Float64[]
  for N in 1:5
    vv = floquet_expansion(H, VanVleck(), N)
    heff = sort!([
      fold(z, wd) for z in eigvals(tomatrix(effective_generator(vv), D, substitutions))
    ])
    push!(errs, maximum(minimum(abs(e - q) for q in Qs) for e in heff))
  end

  @test issorted(errs; rev=true)
  @test errs[end] < 1.0e-7
  @test errs[1] / errs[end] > 1.0e5
end

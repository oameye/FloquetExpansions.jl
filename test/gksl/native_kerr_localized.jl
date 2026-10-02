using Test
using LinearAlgebra

const LocalizedExact = Complex{Rational{BigInt}}

struct LocalizedAffine{T}
  constant::T
  linear::T
end

struct LocalizedFraction{T}
  numerator::LocalizedAffine{T}
  denominator::LocalizedAffine{T}
end

struct LocalizedAffine2{T}
  constant::T
  left::T
  right::T
end

struct LocalizedFraction2{T}
  numerator::LocalizedAffine2{T}
  denominator::LocalizedAffine2{T}
end

localized_value(p::LocalizedAffine, n) = p.constant + p.linear * n
function localized_value(r::LocalizedFraction, n)
  return localized_value(r.numerator, n) / localized_value(r.denominator, n)
end
function localized_value(p::LocalizedAffine2, nleft, nright)
  return p.constant + p.left * nleft + p.right * nright
end
function localized_value(r::LocalizedFraction2, nleft, nright)
  return localized_value(r.numerator, nleft, nright) /
         localized_value(r.denominator, nleft, nright)
end

function localized_kerr_nojump(n, Δ, χ, κ)
  energy = Δ * n + χ * n^2
  decay = κ * n / 2
  return LocalizedExact(-decay, -energy)
end

function localized_kerr_bohr(charge, Δ, χ, κ)
  d0 = localized_kerr_nojump(charge, Δ, χ, κ) - localized_kerr_nojump(0, Δ, χ, κ)
  d1 = localized_kerr_nojump(1 + charge, Δ, χ, κ) - localized_kerr_nojump(1, Δ, χ, κ)
  return LocalizedAffine(d0, d1 - d0)
end

function localized_kerr_bohr2(qleft, qright, Δ, χ, κ)
  function value(nleft, nright)
    left =
      localized_kerr_nojump(nleft + qleft, Δ, χ, κ) - localized_kerr_nojump(nleft, Δ, χ, κ)
    right = conj(
      localized_kerr_nojump(nright + qright, Δ, χ, κ) -
      localized_kerr_nojump(nright, Δ, χ, κ),
    )
    return left + right
  end
  d00 = value(0, 0)
  return LocalizedAffine2(d00, value(1, 0) - d00, value(0, 1) - d00)
end

function localized_creation(d)
  ad = zeros(ComplexF64, d, d)
  for n in 0:(d - 2)
    ad[n + 2, n + 1] = sqrt(n + 1)
  end
  return ad
end

function localized_nojump_superoperator(d, Δ, χ, κ)
  values = ComplexF64[]
  for nright in 0:(d - 1), nleft in 0:(d - 1)
    push!(
      values,
      localized_kerr_nojump(nleft, Δ, χ, κ) + conj(localized_kerr_nojump(nright, Δ, χ, κ)),
    )
  end
  return Diagonal(values)
end

function localized_shift_superoperator(d, qleft, qright, coefficient)
  result = zeros(ComplexF64, d^2, d^2)
  for nright in 0:(d - 1), nleft in 0:(d - 1)
    outleft = nleft + qleft
    outright = nright + qright
    0 <= outleft < d || continue
    0 <= outright < d || continue
    input = nleft + 1 + d * nright
    output = outleft + 1 + d * outright
    result[output, input] = coefficient(nleft, nright)
  end
  return result
end

@testset "exact localized Kerr one-slot charged solve" begin
  Δ = big(1) // big(2)
  χ = big(3) // big(10)
  κ = big(4) // big(5)
  charge = 1

  bohr = localized_kerr_bohr(charge, Δ, χ, κ)
  @test bohr.constant == LocalizedExact(-2 // 5, -4 // 5)
  @test bohr.linear == LocalizedExact(0 // 1, -3 // 5)
  @test !iszero(bohr.linear)

  target = LocalizedAffine(one(LocalizedExact), zero(LocalizedExact))
  solution = LocalizedFraction(target, bohr)

  # A nonconstant affine polynomial cannot divide a nonzero constant in C[n]. The localized
  # coefficient 1/d_1(n), however, solves the charged static equation exactly.
  @test !iszero(target.constant)
  @test iszero(target.linear)
  @test solution.numerator == target
  @test solution.denominator == bohr
  for n in 0:12
    @test localized_value(bohr, n) * localized_value(solution, n) == one(LocalizedExact)
  end

  root = -bohr.constant / bohr.linear
  @test root == LocalizedExact(-4 // 3, 2 // 3)
  @test !iszero(imag(root))
  @test all(!iszero(localized_value(bohr, n)) for n in 0:100)

  # Dense finite-Fock realization of [Qbar, a† r(n)] = a†. This is the same charged static
  # inverse seen by the #370 finite-dimensional oracle, evaluated before any cutoff transition.
  for d in (3, 5, 8)
    occupations = collect(0:(d - 1))
    qdiag = ComplexF64[localized_kerr_nojump(n, Δ, χ, κ) for n in occupations]
    rdiag = ComplexF64[localized_value(solution, n) for n in occupations]
    Qbar = Diagonal(qdiag)
    ad = localized_creation(d)
    X = ad * Diagonal(rdiag)
    @test Qbar * X - X * Qbar ≈ ad atol = 2e-14 rtol = 2e-14
  end
end

@testset "exact localized Kerr two-leg multiplier" begin
  Δ = big(1) // big(2)
  χ = big(3) // big(10)
  κ = big(4) // big(5)
  qleft = -3
  qright = -2

  bohr = localized_kerr_bohr2(qleft, qright, Δ, χ, κ)
  @test bohr.constant == LocalizedExact(2 // 1, -1 // 1)
  @test bohr.left == LocalizedExact(0 // 1, 9 // 5)
  @test bohr.right == LocalizedExact(0 // 1, -6 // 5)

  target = LocalizedAffine2(one(LocalizedExact), zero(LocalizedExact), zero(LocalizedExact))
  solution = LocalizedFraction2(target, bohr)
  @test !iszero(bohr.left)
  @test !iszero(bohr.right)
  @test solution.numerator == target
  @test solution.denominator == bohr

  for nleft in 3:9, nright in 2:8
    νleft = nleft - 3
    νright = nright - 2
    expected = LocalizedExact(5 * κ / 2, Δ + χ * (6 * νleft - 4 * νright + 5))
    @test localized_value(bohr, nleft, nright) == expected
    @test localized_value(bohr, nleft, nright) * localized_value(solution, nleft, nright) ==
      one(LocalizedExact)
  end

  # Re d_{-3,-2}=5κ/2=2 on the whole physical lattice, so the affine zero set is a complex
  # algebraic certificate only; it never becomes a physical Fock resonance for κ>0.
  @test real(bohr.constant) == 2
  @test iszero(real(bohr.left))
  @test iszero(real(bohr.right))
  @test all(!iszero(localized_value(bohr, nleft, nright)) for nleft in 3:30, nright in 2:30)

  for d in (6, 8, 10)
    nojump = localized_nojump_superoperator(d, Δ, χ, κ)
    target_map = localized_shift_superoperator(d, qleft, qright, (_, _) -> 1.0 + 0.0im)
    gauge_map = localized_shift_superoperator(
      d,
      qleft,
      qright,
      (nleft, nright) -> ComplexF64(localized_value(solution, nleft, nright)),
    )
    @test nojump * gauge_map - gauge_map * nojump ≈ target_map atol = 5e-13 rtol = 5e-13
  end
end

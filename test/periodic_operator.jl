using Test
using FloquetExpansions
using LinearAlgebra: ishermitian
import SecondQuantizedAlgebra as SQA
using Symbolics: Differential, expand_derivatives

h = FockSpace(:cavity)
a = Destroy(h, :a)
@variables w::Real t::Real

# `PeriodicOperator` only claims to store `X(t) = ∑_l X_l exp(-i l ω_d t)`. Every test below is
# against that time-dependent operator, reconstructed here and handed to Symbolics' own calculus,
# rather than against the harmonic bookkeeping restated. Each has been checked to FAIL under a
# corrupted implementation (index `p-q` for `p+q`, and a flipped sign in each of d/dt and its
# inverse); a test that cannot fail is not kept.
evaluate(X::PeriodicOperator, t) = X(t)

function ddt(q::SQA.QAdd, t)
  D = Differential(t)
  out = zero(SQA.QAdd)
  for (term, coeff) in q
    mono = isempty(term.ops) ? one(SQA.QAdd) : prod(term.ops)
    out = out + expand_derivatives(D(SQA.to_num(coeff))) * mono
  end
  return out
end

vanishes(q::SQA.QAdd) = iszero(SQA.simplify(q))
vanishes(X::PeriodicOperator) = all(vanishes(X[l]) for l in keys(X))

# A drive with a dense, asymmetric harmonic support, so cross terms genuinely collide in buckets.
function drive()
  return PeriodicOperator(
    Dict(
      0 => a' * a,
      1 => a,
      -1 => a',
      2 => a' * a',
      -2 => a * a,
      3 => a' * a * a,
      -3 => a' * a' * a,
    ),
    w,
  )
end

@testset "commutator is the commutator of the time-dependent operators" begin
  # The convolution ∑_p [K_p, X_{l-p}] is only meaningful if it reproduces [K(t), X(t)].
  # Catches a wrong bucket index and contributions overwritten instead of summed.
  K = drive()
  X = PeriodicOperator(Dict(1 => a, -1 => a', 2 => a' * a'), w)
  @test vanishes(
    evaluate(SQA.commutator(K, X), t) - SQA.commutator(evaluate(K, t), evaluate(X, t))
  )

  # Jacobi, which no amount of per-bucket bookkeeping gets right by accident.
  Y = PeriodicOperator(Dict(1 => a' * a, -2 => a), w)
  jac =
    SQA.commutator(SQA.commutator(K, X), Y) +
    SQA.commutator(SQA.commutator(X, Y), K) +
    SQA.commutator(SQA.commutator(Y, K), X)
  @test vanishes(jac)
end

@testset "derivative is d/dt, in the package's Fourier convention" begin
  # ω_d is factored out of the recursion, so `derivative` is d/dt divided by ω_d.
  X = drive()
  @test vanishes(ddt(evaluate(X, t), t) - w * evaluate(derivative(X), t))
end

@testset "antiderivative inverts d/dt" begin
  X = PeriodicOperator(
    Dict(1 => a, -1 => a', 3 => a' * a, -3 => a' * a, 7 => a' * a' * a), w
  )
  K = antiderivative(X, VanVleck())
  @test vanishes(ddt(evaluate(K, t), t) - w * evaluate(X, t))

  # Inverse in both compositions, which pins the ω_d bookkeeping to be consistent.
  @test vanishes(derivative(K) - X)
  @test vanishes(antiderivative(derivative(X), VanVleck()) - X)
end

@testset "the van Vleck gauge is what fixes the integration constant" begin
  X = drive() - time_average(drive())
  K = antiderivative(X, VanVleck())
  @test vanishes(time_average(K))     # ⟨K⟩ = 0 is the gauge, not an accident of storage

  # The domain is the image of d/dt. A DC harmonic reaching here means a caller forgot the
  # mean subtraction, which is a real bug that once shipped in the draft recursion.
  @test_throws ArgumentError antiderivative(drive(), VanVleck())
end

@testset "Hermiticity survives every operation" begin
  # DESIGN §7 halves the commutator count by filling l < 0 from l >= 0 by adjoint. That is
  # only sound because every triangle node is Hermitian, which needs these four closures.
  K = drive()
  X = PeriodicOperator(Dict(1 => a, -1 => a'), w)
  @test ishermitian(K)
  @test ishermitian(X)
  @test ishermitian(derivative(K))
  @test ishermitian(antiderivative(K - time_average(K), VanVleck()))
  @test ishermitian(im * SQA.commutator(K, X))          # i·ad maps Hermitian to Hermitian
  @test !ishermitian(SQA.commutator(K, X))              # without the i it is anti-Hermitian
end

@testset "harmonic support adds under commutator" begin
  # Underwrites the |l| <= (n+1)M bound the truncation-free claim rests on.
  K = PeriodicOperator(Dict(2 => a, -2 => a'), w)
  X = PeriodicOperator(Dict(3 => a' * a', -3 => a * a), w)
  s = support(SQA.commutator(K, X))
  @test first(s) >= first(support(K)) + first(support(X))
  @test last(s) <= last(support(K)) + last(support(X))
end

@testset "weights stay exact rationals" begin
  # D5: the spec's prefactors (17/4, 3/20, 9/40) and the `iszero(simplify(a-b))` comparison rule
  # both break if a weight silently becomes a float. l = 49 is the smallest harmonic where
  # Float64 fails to round-trip (49 * (1/49) == 0.9999999999999999), so this separates an exact
  # weight from a float one instead of asserting how the coefficient happens to print: exact
  # gives 0, float leaves -1.11e-16.
  K = antiderivative(PeriodicOperator(Dict(49 => a), w), VanVleck())
  @test vanishes(49 * K[49] - im * (1 * a))
end

@testset "inference" begin
  # CLAUDE.md gate. Verified achievable engine-wide: QAdd is concrete and the coefficient
  # Union lives a level below anything these return.
  X = drive()
  Y = PeriodicOperator(Dict(1 => a), w)
  @test @inferred(SQA.commutator(X, Y)) isa PeriodicOperator
  @test @inferred(derivative(X)) isa PeriodicOperator
  @test @inferred(antiderivative(Y, VanVleck())) isa PeriodicOperator
  @test @inferred(time_average(X)) isa SQA.QAdd
  @test @inferred(X + Y) isa PeriodicOperator
  @test @inferred(X - Y) isa PeriodicOperator
  @test @inferred(2 * X) isa PeriodicOperator
  @test @inferred(adjoint(X)) isa PeriodicOperator
  @test @inferred(X[1]) isa SQA.QAdd
end

@testset "the drive frequency is part of the periodic operator" begin
  X = PeriodicOperator(Dict(1 => a, -1 => a'), w)
  @test zero(X) + X == X
  @test_throws MethodError PeriodicOperator(1 => a, -1 => a', w)
  @test_throws ArgumentError X + PeriodicOperator(Dict(1 => a, -1 => a'), 2w)
  @test_throws ArgumentError SQA.commutator(X, PeriodicOperator(Dict(1 => a'), 2w))
end

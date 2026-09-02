using Test
using FloquetExpansions
using LinearAlgebra: ishermitian
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra
using Symbolics: Differential, expand_derivatives

include(joinpath(@__DIR__, "helpers", "shared.jl"))

h = FockSpace(:cavity)
a = Destroy(h, :a)
@variables w::Real t::Real

evaluate(X::PeriodicGenerator, t) = X(t)

function ddt(q::SQA.QAdd, t)
  D = Differential(t)
  out = zero(SQA.QAdd)
  for (term, coeff) in q
    mono = isempty(term.ops) ? one(SQA.QAdd) : prod(term.ops)
    out = out + expand_derivatives(D(SQA.to_num(coeff))) * mono
  end
  return out
end

function drive()
  return PeriodicGenerator(
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
  K = drive()
  X = PeriodicGenerator(Dict(1 => a, -1 => a', 2 => a' * a'), w)
  @test vanishes(
    evaluate(SQA.commutator(K, X), t) - SQA.commutator(evaluate(K, t), evaluate(X, t))
  )

  Y = PeriodicGenerator(Dict(1 => a' * a, -2 => a), w)
  jac =
    SQA.commutator(SQA.commutator(K, X), Y) +
    SQA.commutator(SQA.commutator(X, Y), K) +
    SQA.commutator(SQA.commutator(Y, K), X)
  @test vanishes(jac)
end

@testset "derivative is d/dt, in the package's Fourier convention" begin
  X = drive()
  @test vanishes(ddt(evaluate(X, t), t) - w * evaluate(derivative(X), t))
end

@testset "antiderivative inverts d/dt" begin
  X = PeriodicGenerator(
    Dict(1 => a, -1 => a', 3 => a' * a, -3 => a' * a, 7 => a' * a' * a), w
  )
  K = antiderivative(X, VanVleck())
  @test vanishes(ddt(evaluate(K, t), t) - w * evaluate(X, t))

  @test vanishes(derivative(K) - X)
  @test vanishes(antiderivative(derivative(X), VanVleck()) - X)
end

@testset "the van Vleck gauge is what fixes the integration constant" begin
  X = drive() - time_average(drive())
  K = antiderivative(X, VanVleck())
  @test vanishes(time_average(K))

  @test_throws ArgumentError antiderivative(drive(), VanVleck())
end

@testset "Hermiticity survives every operation" begin
  K = drive()
  X = PeriodicGenerator(Dict(1 => a, -1 => a'), w)
  @test ishermitian(K)
  @test ishermitian(X)
  @test ishermitian(derivative(K))
  @test ishermitian(antiderivative(K - time_average(K), VanVleck()))
  @test ishermitian(im * SQA.commutator(K, X))
  @test !ishermitian(SQA.commutator(K, X))
end

@testset "harmonic support adds under commutator" begin
  K = PeriodicGenerator(Dict(2 => a, -2 => a'), w)
  X = PeriodicGenerator(Dict(3 => a' * a', -3 => a * a), w)
  s = support(SQA.commutator(K, X))
  @test first(s) >= first(support(K)) + first(support(X))
  @test last(s) <= last(support(K)) + last(support(X))
end

@testset "weights stay exact rationals" begin
  K = antiderivative(PeriodicGenerator(Dict(49 => a), w), VanVleck())
  @test vanishes(49 * K[49] - im * (1 * a))
end

@testset "inference" begin
  X = drive()
  Y = PeriodicGenerator(Dict(1 => a), w)
  @test @inferred(SQA.commutator(X, Y)) isa PeriodicGenerator
  @test @inferred(derivative(X)) isa PeriodicGenerator
  @test @inferred(antiderivative(Y, VanVleck())) isa PeriodicGenerator
  @test @inferred(time_average(X)) isa SQA.QAdd
  @test @inferred(X + Y) isa PeriodicGenerator
  @test @inferred(X - Y) isa PeriodicGenerator
  @test @inferred(2 * X) isa PeriodicGenerator
  @test @inferred(adjoint(X)) isa PeriodicGenerator
  @test @inferred(X[1]) isa SQA.QAdd
  @test adjoint(zero(X)) == zero(X)
end

@testset "the drive frequency is part of the periodic generator" begin
  X = PeriodicGenerator(Dict(1 => a, -1 => a'), w)
  @test zero(X) + X == X
  @test_throws MethodError PeriodicGenerator(1 => a, -1 => a', w)
  @test_throws ArgumentError X + PeriodicGenerator(Dict(1 => a, -1 => a'), 2w)
  @test_throws ArgumentError SQA.commutator(X, PeriodicGenerator(Dict(1 => a'), 2w))
  @test_throws MethodError X(0.4)
end

@testset "constructors normalize zero harmonics" begin
  X = PeriodicGenerator(Dict(0 => zero(SQA.QAdd), 1 => a), w)
  @test collect(keys(X)) == [1]
  @test iszero(X[0])
  @test !iszero(X)

  Y = PeriodicGenerator(Dict(0 => zero(SQA.QAdd)), w, zero(SQA.QAdd))
  @test iszero(Y)
  @test Y[0] == zero(SQA.QAdd)

  Z = PeriodicGenerator(Dict{Int,SQA.QAdd}(), w)
  @test iszero(Z)
  @test Z[0] == zero(SQA.QAdd)
  @test_throws ArgumentError PeriodicGenerator(Dict{Int,Any}(), w)
end

using Test
using FloquetExpansions
using LinearAlgebra: ishermitian
import SecondQuantizedAlgebra as SQA

h = FockSpace(:cavity)
a = Destroy(h, :a)
@variables w::Real t::Real g::Real ξ::Real φ::Real t₀::Real

# `harmonics` consumes `exponential_form(H)`, since Base's `exp(::Complex)` expands a written
# `exp(im*w*t)` into cos+i·sin before SQA ever sees it. So the round trip is against that
# normalized form, not against however the caller happened to spell the input.
roundtrips(H) = iszero(SQA.simplify(harmonics(H, w, t)(t) - SQA.exponential_form(H)))

@testset "round trip on a drive that hits every trap at once" begin
  # Merged coefficients, cos normalization, exact rationals, a DC term, and harmonics out to
  # ±2 with operator content that differs per harmonic.
  H =
    (3 // 4) * (a' * a) +
    (1 // 3) * (a' * a' * a * a) +
    g * SQA.expim(-w * t) * (a' * a * a) +
    g * SQA.expim(w * t) * (a' * a' * a) +
    ξ * SQA.expim(-2 * w * t) * (a' * a') +
    ξ * SQA.expim(2 * w * t) * (a * a) +
    cos(w * t) * (a + a')
  @test roundtrips(H)

  X = harmonics(H, w, t)
  @test sort!(collect(keys(X))) == [-2, -1, 0, 1, 2]
  @test ishermitian(X)                       # H_{-m} = H_m†, eq:fourierH

  # And back the other way, so neither direction can be quietly lossy.
  @test harmonics(X(t), w, t) == X
end

@testset "the reported Duffing example collects after rotation" begin
  h_duffing = FockSpace(:duffing_cavity)
  @qnumbers a_duffing::Destroy(h_duffing)
  @variables ω₀_duffing::Real ω_duffing::Real α_duffing::Real F_duffing::Real t_duffing::Real

  x_duffing = (a_duffing + a_duffing') / sqrt(2 * ω₀_duffing)
  p_duffing = im * sqrt(ω₀_duffing / 2) * (a_duffing' - a_duffing)
  duffing_hamiltonian(x, p) =
    ω₀_duffing * (p^2 + x^2) // 2 + α_duffing * x^4 // 4 -
    F_duffing * x * exponential_form(cos(ω_duffing * t_duffing))

  duffing = duffing_hamiltonian(x_duffing, p_duffing)
  rotation = Rotation(a_duffing, ω_duffing * t_duffing, t_duffing)
  H = transform(duffing, rotation) |> simplify
  X = harmonics(H, ω_duffing, t_duffing)

  @test sort!(collect(keys(X))) == [-4, -2, 0, 2, 4]
  @test iszero(SQA.simplify(X(t_duffing) - SQA.exponential_form(H)))
end

@testset "coefficients that merge across harmonics are still split" begin
  # THE finding that broke the naive design (DESIGN §3 F1). A QAdd maps monomial -> coefficient,
  # so these two terms share the monomial `a` and are stored as ONE entry whose coefficient is a
  # sum over two different harmonics. Classifying per stored term rather than per additive part
  # silently drops one of them.
  H = 2 * SQA.expim(-w * t) * a + 5 * SQA.expim(2 * w * t) * a
  @test length(H) == 1                       # one stored entry, two harmonics inside it

  X = harmonics(H, w, t)
  @test sort!(collect(keys(X))) == [-2, 1]
  @test iszero(SQA.simplify(X[1] - 2 * a))
  @test iszero(SQA.simplify(X[-2] - 5 * a))
  @test roundtrips(H)
end

@testset "a phase behind a real denominator is still split" begin
  # This enters SQA's raw coefficient tier; phase decomposition must survive that boundary.
  H = ((-0.5 + 0im) * g + (-0.5 + 0im) * g * SQA.expim(-2 * w * t)) / sqrt(2 * ξ) * a
  X = harmonics(H, w, t)
  @test sort!(collect(keys(X))) == [0, 2]
  @test roundtrips(H)
end

@testset "trigonometric input is normalized before parsing" begin
  # A written `exp(im*w*t)` never reaches SQA intact, so cos/sin is the form users actually
  # hand over. Both must land on ±1 with the right weights.
  Xc = harmonics(cos(w * t) * (a + a'), w, t)
  @test sort!(collect(keys(Xc))) == [-1, 1]
  @test iszero(SQA.simplify(Xc[1] - Xc[-1]))          # cos is even
  @test roundtrips(cos(w * t) * (a + a'))

  Xs = harmonics(im * sin(w * t) * (a - a'), w, t)
  @test sort!(collect(keys(Xs))) == [-1, 1]
  @test iszero(SQA.simplify(Xs[1] + Xs[-1]))          # sin is odd
  @test roundtrips(im * sin(w * t) * (a - a'))
end

@testset "a constant phase offset belongs in the coefficient, not the index" begin
  # cos(w*t + φ) is an ordinary drive. The φ must not be read as part of the harmonic index,
  # and must come back on the coefficient as exp(∓iφ).
  H = cos(w * t + φ) * (a + a')
  X = harmonics(H, w, t)
  @test sort!(collect(keys(X))) == [-1, 1]
  @test roundtrips(H)
  @test ishermitian(X)
end

@testset "a frequency-dependent static phase offset is retained" begin
  H = SQA.expim(w * (t + t₀)) * a
  X = harmonics(H, w, t)
  @test collect(keys(X)) == [-1]
  @test iszero(SQA.simplify(X[-1] - SQA.expim(w * t₀) * a))
  @test roundtrips(H)
end

@testset "a phase over a denominator still finds its harmonic" begin
  # A static symbolic denominator belongs to the amplitude, not the harmonic index.
  H = (g / w) * SQA.expim(-w * t) * a + conj(g / w) * SQA.expim(w * t) * a'
  X = harmonics(H, w, t)
  @test sort!(collect(keys(X))) == [-1, 1]
  @test iszero(SQA.simplify(X[1] - (g / w) * a))
  @test roundtrips(H)
  @test ishermitian(X)
end

@testset "non-periodic and fractional phases are rejected" begin
  # Silently accepting either would produce a wrong expansion rather than an error.
  @test_throws ArgumentError harmonics(SQA.expim(w * t^2) * a, w, t)
  @test_throws ArgumentError harmonics(SQA.expim((1 // 2) * w * t) * a, w, t)
end

@testset "exact coefficients survive the collector" begin
  # D5 again: the parser must not be where 1//3 turns into 0.333…, or the engine's exact
  # comparison against the spec's prefactors is lost before it starts.
  X = harmonics((1 // 3) * (a' * a) + (1 // 7) * SQA.expim(-w * t) * a, w, t)
  @test iszero(SQA.simplify(3 * X[0] - 1 * (a' * a)))
  @test iszero(SQA.simplify(7 * X[1] - 1 * a))
end

@testset "a static operator is a single DC harmonic" begin
  X = harmonics(1 * (a' * a), w, t)
  @test collect(keys(X)) == [0]
  @test iszero(SQA.simplify(time_average(X) - 1 * (a' * a)))
  @test iszero(derivative(X))
end

@testset "PeriodicOperator accepts symbolic drives directly" begin
  H = a * SQA.expim(w * t) + a' * SQA.expim(-w * t)
  X = PeriodicOperator(H, w)
  @test X == harmonics(H, w, t)
  @test @inferred(PeriodicOperator(H, w)) isa PeriodicOperator
  @test @inferred(PeriodicOperator(H, w, t)) isa PeriodicOperator
  @test iszero(SQA.simplify(X(t) - SQA.exponential_form(H)))

  static = PeriodicOperator(a' * a, w)
  @test collect(keys(static)) == [0]
  @test iszero(SQA.simplify(static[0] - a' * a))
end

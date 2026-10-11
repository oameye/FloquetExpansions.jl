using Test
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions
const Q = Rational{BigInt}

function synthetic(x::Vector{Q})
  x1, x2, x3 = x
  return Dict{Symbol,Q}(
    :mixed => 2 * x1^2 * x2 - 3 // 7 * x2^3 + x1 * x3^2 + 5,
    :constant => 5 // 3,
    :homogeneous => x1^5 * x2 * x3,
    :imaginary => -x2 * x3,
  )
end

const SYNTHETIC_KEYS = collect(keys(synthetic(ones(Q, 3))))

@testset "polynomial reconstruction recovers sparse polynomials in each sign region" begin
  for signs in ([1, 1, 1], [-1, 1, -1])
    functions, verified = FE.reconstruct_polynomials(synthetic, SYNTHETIC_KEYS, signs, 3)
    @test length(verified) == 3
    @test all(point -> sign.(point) == signs, verified)
    @test length(functions[:mixed]) == 4
    @test length(functions[:homogeneous]) == 1
    for magnitudes in
        (Q[3 // 4, 5 // 2, 7 // 3], Q[2 // 9, 11 // 5, 1 // 6], Q[40, 1 // 30, 9])
      point = signs .* magnitudes
      expected = synthetic(point)
      for key in SYNTHETIC_KEYS
        @test FE.polynomial_value(functions[key], point) == expected[key]
      end
    end
  end
end

@testset "polynomial reconstruction is bound to the declared sign region" begin
  odd(x::Vector{Q}) = Dict{Symbol,Q}(:f => x[1] < 0 ? -x[1]^3 : x[1]^3)
  positive, _ = FE.reconstruct_polynomials(odd, [:f], [1], 2)
  negative, _ = FE.reconstruct_polynomials(odd, [:f], [-1], 2)
  @test FE.polynomial_value(positive[:f], Q[2]) == 8
  @test FE.polynomial_value(negative[:f], Q[-2]) == 8
  @test FE.polynomial_value(positive[:f], Q[-2]) != odd(Q[-2])[:f]
end

@testset "polynomial reconstruction rejects a branch beyond the magnitudes it fitted at" begin
  branch(limit) = x::Vector{Q} -> Dict{Symbol,Q}(:f => x[1] < limit ? x[1]^2 : x[1]^2 + 1)
  @test_throws ArgumentError FE.reconstruct_polynomials(branch(12), [:f], [1], 8)
  functions, _ = FE.reconstruct_polynomials(branch(10^6), [:f], [1], 8)
  @test FE.polynomial_value(functions[:f], Q[3]) == 9
end

@testset "polynomial reconstruction reports what it cannot fit" begin
  dense(x::Vector{Q}) = Dict{Symbol,Q}(:f => sum(x)^12)
  @test_throws ArgumentError FE.reconstruct_polynomials(dense, [:f], ones(Int, 4), 1)
  steep(x::Vector{Q}) = Dict{Symbol,Q}(:f => x[1]^45)
  @test_throws ArgumentError FE.reconstruct_polynomials(steep, [:f], [1], 1)
  ratio(x::Vector{Q}) = Dict{Symbol,Q}(:f => x[1] / (1 + x[1]^2))
  @test_throws ArgumentError FE.reconstruct_polynomials(ratio, [:f], [1], 1)
end

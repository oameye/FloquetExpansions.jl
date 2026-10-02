using Test
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions

function synthetic_outputs(x::Vector{Rational{BigInt}})
  x1, x2, x3 = x
  return Dict{Symbol,FE.ExactField}(
    :homogeneous => 2 * x1^2 * x2 - 3 // 7 * x2^3 + im * x1 * x3^2,
    :shifted => (x1 + 2 * x2^2) / (1 + x1 * x3),
    :ratio => x1 / x2,
    :constant => 5 // 3 + 0im,
  )
end

@testset "exact rational reconstruction recovers sparse rational functions" begin
  for signs in ([1, 1, 1], [-1, 1, 1])
    sampler = FE.ExactSampler(
      synthetic_outputs, collect(keys(synthetic_outputs(ones(Rational{BigInt}, 3)))), signs
    )
    fitted = FE.reconstruct_rational_functions(sampler)
    @test length(fitted[:homogeneous].numerator) == 3
    @test length(fitted[:homogeneous].denominator) == 1
    @test length(fitted[:shifted].denominator) == 2
    @test length(fitted[:constant].numerator) == 1
    for values in ([3 // 4, 5 // 2, 7 // 3], [2 // 9, 11 // 5, 1 // 6])
      point = signs .* Rational{BigInt}.(values)
      expected = synthetic_outputs(point)
      for (key, f) in fitted
        @test FE.rational_value(f, FE.ExactField.(point)) == expected[key]
      end
    end
  end
end

@testset "exact rational reconstruction rejects a point that breaks the structure" begin
  function piecewise(x::Vector{Rational{BigInt}})
    return Dict{Symbol,FE.ExactField}(:f => x[1] < 1 ? x[1]^2 : x[1]^2 + x[1]^5)
  end
  sampler = FE.ExactSampler(piecewise, [:f], [1])
  @test_throws ArgumentError FE.reconstruct_rational_functions(sampler)
end

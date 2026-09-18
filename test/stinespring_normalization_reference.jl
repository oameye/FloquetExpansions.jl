using Test
using LinearAlgebra: I, inv

include(joinpath(@__DIR__, "helpers", "stinespring_normalization_reference.jl"))

const StinespringExact = Rational{Int}

function exact_identity(n::Int)
  return Matrix{StinespringExact}(I, n, n)
end

@testset "physical Kraus keys control left/right pasting" begin
  I2 = exact_identity(2)
  A = StinespringExact[1 1; 0 1]
  B = StinespringExact[0 1; 1 0]
  C = StinespringExact[1 0; 1 -1]

  histories = [
    StinespringHistory(0, :vacuum, (:vacuum,), [I2]),
    StinespringHistory(1, :shared_output, (:history_a,), [A]),
    StinespringHistory(1, :shared_output, (:history_b,), [B]),
    StinespringHistory(1, :orthogonal_output, (:history_c,), [C]),
  ]
  rows = coalesce_stinespring_histories(histories)
  metric = stinespring_metric_coefficients(rows, 1, I2)

  expected = adjoint(A + B) * (A + B) + adjoint(C) * C
  naive = adjoint(A) * A + adjoint(B) * B + adjoint(C) * C

  @test metric[1] == I2
  @test metric[2] == expected
  @test expected != naive
  @test length(rows) == 3
end

@testset "finite amplitude triangle has only a beyond-order TP defect" begin
  I2 = exact_identity(2)
  U = StinespringExact[0 1; 1 0]

  for generator_order in 0:5
    rows = stinespring_cayley_rows(generator_order, U)
    metric_order = generator_order + 3
    metric = stinespring_metric_coefficients(rows, metric_order, I2)

    @test metric[1] == I2
    for r in 1:(generator_order + 1)
      @test iszero(metric[r + 1])
    end
    @test !iszero(metric[generator_order + 3])

    inverse_sqrt = stinespring_inverse_sqrt_series(metric, metric_order, I2)
    for r in 1:(generator_order + 1)
      @test iszero(inverse_sqrt[r + 1])
    end

    normalized = stinespring_normalize_rows_series(
      rows, inverse_sqrt, generator_order + 1, I2
    )
    for (row, normalized_row) in zip(rows, normalized)
      for n in 0:(length(row.coefficients) - 1)
        stinespring_triangle_contains(generator_order, row.output_number, n) || continue
        @test normalized_row.coefficients[n + 1] == row.coefficients[n + 1]
      end
    end
  end
end

@testset "noncommutative inverse-square-root recurrence selects the principal branch" begin
  I2 = exact_identity(2)
  Z2 = zero(I2)
  A = StinespringExact[1 1; 1 0]
  B = StinespringExact[0 1; 1 2]
  @test A * B != B * A

  positive_factor = [I2, A, B, Z2, Z2]
  metric = stinespring_series_product(positive_factor, positive_factor, 4, Z2)
  inverse_sqrt = stinespring_inverse_sqrt_series(metric, 4, I2)
  expected_inverse = stinespring_inverse_series(positive_factor, 4, I2)

  @test inverse_sqrt == expected_inverse

  left_product = stinespring_series_product(inverse_sqrt, metric, 4, Z2)
  normalized_metric = stinespring_series_product(left_product, inverse_sqrt, 4, Z2)
  @test normalized_metric[1] == I2
  @test all(iszero, normalized_metric[2:end])
  @test all(coefficient == adjoint(coefficient) for coefficient in inverse_sqrt)
end

@testset "finite right-normalization is exactly isometric" begin
  I2 = exact_identity(2)
  U = StinespringExact[0 1; 1 0]
  K0 = (3 // 5) * I2
  K1 = (4 // 5) * U
  @test adjoint(K0) * K0 + adjoint(K1) * K1 == I2

  positive_factor = StinespringExact[1 1//5; 1//5 6//5]
  @test positive_factor == adjoint(positive_factor)
  right_normalizer = inv(positive_factor)

  V0 = K0 * positive_factor
  V1 = K1 * positive_factor
  completeness = adjoint(V0) * V0 + adjoint(V1) * V1
  @test completeness == positive_factor * positive_factor
  @test right_normalizer * completeness * right_normalizer == I2

  normalized0 = V0 * right_normalizer
  normalized1 = V1 * right_normalizer
  @test normalized0 == K0
  @test normalized1 == K1
  @test adjoint(normalized0) * normalized0 + adjoint(normalized1) * normalized1 == I2
end

using Test
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions

@testset "graded channel Gram coefficients are exact Cauchy products" begin
  T = Complex{Rational{Int}}
  r0 = T[1, 0, 2 // 3]
  r1 = T[0, 1 // 2 + im, -1]
  r2 = T[3, -im // 4, 0]
  u0 = T[1 // 5, 2im, 1]
  channels = [FE.GradedChannel{T}(0, [r0, r1, r2]), FE.GradedChannel{T}(1, [u0])]

  outer(x, y) = x * y'
  hermitian(X) = (X + X') / 2
  @test @inferred(FE.gram_coefficient(channels, 0, 3)) == outer(r0, r0)
  # The order-one rate of the channel born at order one is its leading square.
  @test FE.gram_coefficient(channels, 1, 3) == hermitian(2 * outer(r1, r0)) + outer(u0, u0)
  @test FE.gram_coefficient(channels, 2, 3) == hermitian(2 * outer(r2, r0)) + outer(r1, r1)
  @test eltype(FE.gram_coefficient(channels, 2, 3)) == T

  # K_n^< excludes every product containing a leading amplitude.
  @test @inferred(FE.known_gram(channels, 2, 3)) == outer(r1, r1)
  @test FE.known_gram(channels[1:1], 3, 3) == hermitian(2 * outer(r2, r1))
  # The order-one channel still lacks its first correction, so K_3^< is not yet defined.
  @test_throws ArgumentError FE.known_gram(channels, 3, 3)
  @test iszero(FE.known_gram(channels, 1, 3))
  @test_throws ArgumentError FE.known_gram(channels, 5, 3)
end

@testset "graded channels activate after birth and stay prefix complete" begin
  T = ComplexF64
  channels = [FE.GradedChannel{T}(0, [T[1, 0]])]
  indices, active = @inferred FE.active_channels(channels, 1, 2)
  @test indices == [1]
  @test active == reshape(T[1, 0], :, 1)

  FE.store_births!(channels, reshape(T[0, 2], :, 1), 1)
  @test [channel.onset for channel in channels] == [0, 1]
  # A newborn is not active at its own birth order.
  @test first(FE.active_channels(channels, 1, 2)) == [1]
  @test first(FE.active_channels(channels, 2, 2)) == [1, 2]

  FE.store_corrections!(channels, [1], reshape(T[0, 1], :, 1), 1)
  @test channels[1].coefficients == [T[1, 0], T[0, 1]]
  @test_throws ArgumentError FE.store_corrections!(channels, [1], reshape(T[0, 1], :, 1), 1)
  @test_throws DimensionMismatch FE.store_corrections!(channels, [1, 2], zeros(T, 2, 1), 2)

  empty_indices, empty_active = FE.active_channels(FE.GradedChannel{T}[], 3, 4)
  @test isempty(empty_indices)
  @test size(empty_active) == (4, 0)
end

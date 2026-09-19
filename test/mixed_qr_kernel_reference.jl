using Test

include(joinpath(@__DIR__, "helpers", "simplex_output_kernel_reference.jl"))
include(joinpath(@__DIR__, "helpers", "mixed_qr_kernel_reference.jl"))

@testset "formal-period simplex recurrence retains folded pieces exactly" begin
  primitive_two = mixed_formal_period_simplex((1, -1))
  @test mixed_period_coefficient(primitive_two, 1) == mixed_kernel_im
  @test mixed_period_coefficient(primitive_two, 2) == 0

  primitive_three = mixed_formal_period_simplex((2, -1, -1))
  @test mixed_period_coefficient(primitive_three, 1) == -(1 // 2)
  @test mixed_period_coefficient(primitive_three, 2) == 0

  reducible = mixed_formal_period_simplex((1, -1, 0))
  @test mixed_period_coefficient(reducible, 1) == 1
  @test mixed_period_coefficient(reducible, 2) == mixed_kernel_im / 2
end

@testset "bare physical one-output wavepackets are period orthogonal" begin
  matched = mixed_qr_gram_polynomial(
    MixedQRVertex[mixed_qr_jump(:a, 2)], MixedQRVertex[mixed_qr_jump(:a, 2)]
  )
  @test mixed_period_coefficient(matched, 1) == 1

  sideband_orthogonal = mixed_qr_gram_polynomial(
    MixedQRVertex[mixed_qr_jump(:a, 2)], MixedQRVertex[mixed_qr_jump(:a, 1)]
  )
  @test all(iszero, sideband_orthogonal)

  channel_orthogonal = mixed_qr_gram_polynomial(
    MixedQRVertex[mixed_qr_jump(:a, 2)], MixedQRVertex[mixed_qr_jump(:b, 2)]
  )
  @test all(iszero, channel_orthogonal)
end

@testset "one drift insertion closes in a two-wavepacket physical subspace" begin
  drift_harmonic = 2
  jump_harmonic = 3
  inverse_frequency = one(MixedKernelExact) / (mixed_kernel_im * drift_harmonic)

  before = mixed_one_drift_jump_kernel(:before, drift_harmonic, jump_harmonic)
  @test length(before) == 2
  @test before[jump_harmonic] == inverse_frequency
  @test before[jump_harmonic + drift_harmonic] == -inverse_frequency

  after = mixed_one_drift_jump_kernel(:after, drift_harmonic, jump_harmonic)
  @test length(after) == 2
  @test after[jump_harmonic] == -inverse_frequency
  @test after[jump_harmonic + drift_harmonic] == inverse_frequency
  @test after == Dict(harmonic => -coefficient for (harmonic, coefficient) in before)
end

@testset "mixed Q/R Gram sum gives exact dressed-jump norm and interference" begin
  drift_harmonic = 2
  jump_harmonic = 3
  expected_norm = MixedKernelExact(2 // drift_harmonic^2, 0 // 1)

  before = MixedQRVertex[mixed_qr_drift(drift_harmonic), mixed_qr_jump(:a, jump_harmonic)]
  after = MixedQRVertex[mixed_qr_jump(:a, jump_harmonic), mixed_qr_drift(drift_harmonic)]

  before_poset = mixed_qr_poset(before, before)
  after_poset = mixed_qr_poset(after, after)
  @test length(mixed_qr_linear_extensions(before_poset.predecessors)) == 2
  @test length(mixed_qr_linear_extensions(after_poset.predecessors)) == 2

  before_norm = mixed_qr_gram_polynomial(before, before)
  after_norm = mixed_qr_gram_polynomial(after, after)
  cross = mixed_qr_gram_polynomial(before, after)

  @test mixed_period_coefficient(before_norm, 1) == expected_norm
  @test mixed_period_coefficient(after_norm, 1) == expected_norm
  @test mixed_period_coefficient(cross, 1) == -expected_norm
  @test all(iszero(mixed_period_coefficient(before_norm, degree)) for degree in 2:3)
  @test all(iszero(mixed_period_coefficient(after_norm, degree)) for degree in 2:3)
  @test all(iszero(mixed_period_coefficient(cross, degree)) for degree in 2:3)
end

@testset "asymmetric mixed branch words preserve exact Gram hermiticity" begin
  one_left = MixedQRVertex[mixed_qr_drift(1), mixed_qr_jump(:a, 2), mixed_qr_drift(-2)]
  one_right = MixedQRVertex[mixed_qr_jump(:a, 1), mixed_qr_drift(3)]
  one_forward = mixed_qr_gram_polynomial(one_left, one_right)
  one_reverse = mixed_qr_gram_polynomial(one_right, one_left)

  @test length(
    mixed_qr_linear_extensions(mixed_qr_poset(one_left, one_right).predecessors)
  ) == 2
  @test mixed_period_coefficient(one_forward, 1) == mixed_kernel_im / 6
  @test mixed_period_coefficient(one_reverse, 1) == -mixed_kernel_im / 6
  @test all(
    mixed_period_coefficient(one_forward, degree) ==
    conj(mixed_period_coefficient(one_reverse, degree)) for degree in 0:4
  )

  two_left = MixedQRVertex[mixed_qr_drift(1), mixed_qr_jump(:a, 1), mixed_qr_jump(:b, -1)]
  two_right = MixedQRVertex[mixed_qr_jump(:a, 0), mixed_qr_drift(2), mixed_qr_jump(:b, -2)]
  two_forward = mixed_qr_gram_polynomial(two_left, two_right)
  two_reverse = mixed_qr_gram_polynomial(two_right, two_left)

  @test mixed_period_coefficient(two_forward, 1) == -(3 // 2) * mixed_kernel_im
  @test mixed_period_coefficient(two_reverse, 1) == (3 // 2) * mixed_kernel_im
  @test all(
    mixed_period_coefficient(two_forward, degree) ==
    conj(mixed_period_coefficient(two_reverse, degree)) for degree in 0:4
  )
end

@testset "pure-jump mixed reference reduces to the certified simplex overlap" begin
  left_word = SimplexOutputWord((:a, :a, :a), (2, -1, -1))
  right_word = SimplexOutputWord((:a, :a, :a), (0, 0, 0))
  simplex_overlap = simplex_slow_overlap(
    left_word,
    right_word;
    imaginary=mixed_kernel_im,
    inverse_weight=harmonic -> 1 // harmonic,
  )
  @test simplex_overlap.class == SimplexSlowPrimitive

  left = MixedQRVertex[mixed_qr_jump(:a, 2), mixed_qr_jump(:a, -1), mixed_qr_jump(:a, -1)]
  right = MixedQRVertex[mixed_qr_jump(:a, 0), mixed_qr_jump(:a, 0), mixed_qr_jump(:a, 0)]
  gram = mixed_qr_gram_polynomial(left, right)

  @test mixed_period_coefficient(gram, 1) == simplex_overlap.coefficient
  @test mixed_period_coefficient(gram, 2) == 0
end

@testset "fixed delta order gives a finite mixed physical word set" begin
  @test mixed_qr_word_count(6, 3, 0) == 1
  @test mixed_qr_word_count(6, 3, 4) == 1927
  @test mixed_qr_word_count(6, 3, 6) > mixed_qr_word_count(6, 3, 4)
end

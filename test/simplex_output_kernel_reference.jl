using Test
using LinearAlgebra: I, kron

include(joinpath(@__DIR__, "helpers", "simplex_output_kernel_reference.jl"))

const SimplexExact = Complex{Rational{Int}}
const simplex_im = SimplexExact(0 // 1, 1 // 1)

function simplex_matrix_product(left::Matrix{SimplexExact}, right::Matrix{SimplexExact})
  return left * right
end

function simplex_recycling_super(left::Matrix{SimplexExact}, right::Matrix{SimplexExact})
  return kron(conj.(right), left)
end

function simplex_hamiltonian_super(H::Matrix{SimplexExact})
  identity = Matrix{SimplexExact}(I, 2, 2)
  return -simplex_im * (kron(identity, H) - kron(transpose(H), identity))
end

@testset "finite physical simplex kernels replace output Fourier enumeration" begin
  one_left = SimplexOutputWord((:a,), (2,))
  one_right = SimplexOutputWord((:a,), (2,))
  one_overlap = simplex_slow_overlap(
    one_left, one_right; imaginary=simplex_im, inverse_weight=harmonic -> 1 // harmonic
  )
  @test simplex_output_number(one_left) == 1
  @test simplex_total_sideband(one_left) == 2
  @test one_overlap.class == SimplexSlowPrimitive
  @test one_overlap.coefficient == 1

  outside = simplex_slow_overlap(
    SimplexOutputWord((:a,), (2,)),
    SimplexOutputWord((:a,), (1,));
    imaginary=simplex_im,
    inverse_weight=harmonic -> 1 // harmonic,
  )
  @test outside.class == SimplexOutsideSlowSector

  orthogonal = simplex_slow_overlap(
    SimplexOutputWord((:a,), (2,)),
    SimplexOutputWord((:b,), (2,));
    imaginary=simplex_im,
    inverse_weight=harmonic -> 1 // harmonic,
  )
  @test orthogonal.class == SimplexChannelOrthogonal

  @test simplex_output_word_count(2, 3, 4) == 6^4
  @test simplex_output_triangle_count(2, 3, 3) == sum(6^q for q in 0:4)
end

@testset "simplex overlap is the charged first-return coefficient" begin
  two_left = SimplexOutputWord((:a, :a), (1, -1))
  two_right = SimplexOutputWord((:a, :a), (0, 0))
  two_overlap = simplex_slow_overlap(
    two_left, two_right; imaginary=simplex_im, inverse_weight=harmonic -> 1 // harmonic
  )
  @test two_overlap.class == SimplexSlowPrimitive
  @test two_overlap.coefficient == simplex_im

  two_reverse = simplex_slow_overlap(
    two_right, two_left; imaginary=simplex_im, inverse_weight=harmonic -> 1 // harmonic
  )
  @test two_reverse.class == SimplexSlowPrimitive
  @test two_reverse.coefficient == -simplex_im

  reducible = simplex_slow_overlap(
    SimplexOutputWord((:a, :a, :a), (1, -1, 0)),
    SimplexOutputWord((:a, :a, :a), (0, 0, 0));
    imaginary=simplex_im,
    inverse_weight=harmonic -> 1 // harmonic,
  )
  @test reducible.class == SimplexSlowReducible

  three_overlap = simplex_slow_overlap(
    SimplexOutputWord((:a, :a, :a), (2, -1, -1)),
    SimplexOutputWord((:a, :a, :a), (0, 0, 0));
    imaginary=simplex_im,
    inverse_weight=harmonic -> 1 // harmonic,
  )
  @test three_overlap.class == SimplexSlowPrimitive
  @test three_overlap.coefficient == -(1 // 2)
end

@testset "simplex-kernel Gram pasting reconstructs rotating-jump RR" begin
  sigma_x = SimplexExact[0 1; 1 0]
  sigma_y = SimplexExact[0 -simplex_im; simplex_im 0]
  sigma_z = SimplexExact[1 0; 0 -1]
  sigma_minus = (1 // 2) * (sigma_x - simplex_im * sigma_y)
  sigma_plus = (1 // 2) * (sigma_x + simplex_im * sigma_y)
  amplitudes = Dict(
    0 => sigma_z,
    1 => simplex_im * sigma_minus,
    -1 => -simplex_im * sigma_plus,
  )
  zero_super = zeros(SimplexExact, 4, 4)

  rr1 = simplex_two_event_slow_paste(
    amplitudes;
    channel=:jump,
    product=simplex_matrix_product,
    paste=simplex_recycling_super,
    imaginary=simplex_im,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component=zero_super,
    mismatch_abs=1,
  )
  rr2 = simplex_two_event_slow_paste(
    amplitudes;
    channel=:jump,
    product=simplex_matrix_product,
    paste=simplex_recycling_super,
    imaginary=simplex_im,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component=zero_super,
    mismatch_abs=2,
  )
  rr = simplex_two_event_slow_paste(
    amplitudes;
    channel=:jump,
    product=simplex_matrix_product,
    paste=simplex_recycling_super,
    imaginary=simplex_im,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component=zero_super,
  )

  @test rr1 == simplex_hamiltonian_super(-sigma_z)
  @test rr2 == simplex_hamiltonian_super(-(1 // 4) * sigma_z)
  @test rr == rr1 + rr2
  @test rr == simplex_hamiltonian_super(-(5 // 4) * sigma_z)
end

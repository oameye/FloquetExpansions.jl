using Test
using LinearAlgebra: kron

include(joinpath(@__DIR__, "helpers", "open_leg_reference.jl"))

const ExactComplex = Complex{Rational{Int}}
const exact_im = ExactComplex(0 // 1, 1 // 1)

matrix_product(left::Matrix{ExactComplex}, right::Matrix{ExactComplex}) = left * right
matrix_commutator(left::Matrix{ExactComplex}, right::Matrix{ExactComplex}) =
  left * right - right * left

I2_exact = ExactComplex[1 0; 0 1]
Z2_exact = zeros(ExactComplex, 2, 2)
σx_exact = ExactComplex[0 1; 1 0]
σy_exact = ExactComplex[0 -exact_im; exact_im 0]
σz_exact = ExactComplex[1 0; 0 -1]

@testset "open-leg Bloch-Feshbach and Hori-Deprit reference recurrences" begin
  M0 = ExactComplex[1 2; 0 -1]
  M1 = ExactComplex[0 1 + exact_im; 2 -1]
  Mm1 = ExactComplex[1 - exact_im 0; 1 2]
  M2 = ExactComplex[2 0; 1 exact_im]
  Mm2 = ExactComplex[-1 1; -exact_im 2]

  N0 = ExactComplex[1 0; 2 -2]
  N1 = ExactComplex[0 1; -1 1]
  Nm1 = ExactComplex[2 -1; 0 -1]

  A1 = openleg_periodic(
    Dict(
      (0, 1) => M0,
      (1, 1) => M1,
      (-1, 1) => Mm1,
      (2, 1) => M2,
      (-2, 1) => Mm2,
    ),
    Z2_exact,
  )
  A2 = openleg_periodic(
    Dict((0, 0) => N0, (1, 0) => N1, (-1, 0) => Nm1), Z2_exact
  )

  bloch = openleg_bloch_reference(
    [A1, A2], 4; product=matrix_product, identity_component=I2_exact
  )
  lower = openleg_bloch_reference(
    [A1, A2], 3; product=matrix_product, identity_component=I2_exact
  )
  hori = openleg_hori_deprit_order2(A1, A2; product=matrix_product)

  @test openleg_equal(bloch.effective[1], hori.effective1)
  @test openleg_equal(bloch.effective[2], hori.effective2)
  @test openleg_equal(bloch.wave[1], hori.generator1)

  expected_two_output = zero(Z2_exact)
  for harmonic in (1, 2)
    expected_two_output +=
      (-exact_im * (1 // harmonic)) *
      matrix_commutator(A1[harmonic, 1], A1[-harmonic, 1])
  end
  @test bloch.effective[2][0, 2] == expected_two_output
  @test bloch.effective[2][0, 0] == A2[0, 0]

  for n in 1:4
    allowed = Set(q for q in 0:n if iseven(n - q))
    @test all(grade in allowed for grade in openleg_grades(bloch.effective[n]))
    @test all(grade in allowed for grade in openleg_grades(bloch.wave[n]))
  end

  for n in 1:3
    @test openleg_equal(bloch.effective[n], lower.effective[n])
    @test openleg_equal(bloch.wave[n], lower.wave[n])
  end
end

function recycling_super(left::Matrix{ExactComplex}, right::Matrix{ExactComplex})
  return kron(conj.(right), left)
end

function cross_dissipator_super(
  left::Matrix{ExactComplex}, right::Matrix{ExactComplex}
)
  product = adjoint(right) * left
  recycling = recycling_super(left, right)
  anticommutator = kron(I2_exact, product) + kron(transpose(product), I2_exact)
  return recycling - (1 // 2) * anticommutator
end

function hamiltonian_super(H::Matrix{ExactComplex})
  return -exact_im * (kron(I2_exact, H) - kron(transpose(H), I2_exact))
end

function amplitude_pair_harmonic(amplitudes, harmonic::Int, builder)
  result = zeros(ExactComplex, 4, 4)
  for (left_harmonic, left) in amplitudes
    for (right_harmonic, right) in amplitudes
      left_harmonic - right_harmonic == harmonic || continue
      result += builder(left, right)
    end
  end
  return result
end

function loss_harmonic(amplitudes, harmonic::Int)
  result = zero(I2_exact)
  for (left_harmonic, left) in amplitudes
    for (right_harmonic, right) in amplitudes
      left_harmonic - right_harmonic == harmonic || continue
      result += adjoint(right) * left
    end
  end
  return result
end

@testset "rotating-jump open-leg history oracle" begin
  σminus = (1 // 2) * (σx_exact - exact_im * σy_exact)
  σplus = (1 // 2) * (σx_exact + exact_im * σy_exact)
  amplitudes = Dict(
    0 => σz_exact,
    1 => exact_im * σminus,
    -1 => -exact_im * σplus,
  )

  @test loss_harmonic(amplitudes, 0) == 2 * I2_exact
  for harmonic in (-2, -1, 1, 2)
    @test loss_harmonic(amplitudes, harmonic) == zero(I2_exact)
  end

  dissipative = Dict(
    harmonic => amplitude_pair_harmonic(amplitudes, harmonic, cross_dissipator_super) for
    harmonic in -2:2
  )
  recycling = Dict(
    harmonic => amplitude_pair_harmonic(amplitudes, harmonic, recycling_super) for
    harmonic in -2:2
  )

  for harmonic in (-2, -1, 1, 2)
    @test dissipative[harmonic] == recycling[harmonic]
  end

  averaged_rows =
    cross_dissipator_super(amplitudes[0], amplitudes[0]) +
    cross_dissipator_super(amplitudes[1], amplitudes[1]) +
    cross_dissipator_super(amplitudes[-1], amplitudes[-1])
  @test dissipative[0] == averaged_rows

  rr1 = -exact_im * matrix_commutator(dissipative[1], dissipative[-1])
  rr2 =
    (-exact_im * (1 // 2)) * matrix_commutator(dissipative[2], dissipative[-2])
  rr = rr1 + rr2

  @test rr1 == hamiltonian_super(-σz_exact)
  @test rr2 == hamiltonian_super(-(1 // 4) * σz_exact)
  @test rr == hamiltonian_super(-(5 // 4) * σz_exact)

  # The nonzero m=±2 density harmonics are generated solely by pasting opposite
  # one-branch m=±1 sidebands; the one-branch amplitude itself has no ±2 component.
  @test !all(iszero, recycling[2])
  @test !all(iszero, recycling[-2])
end

@testset "connected two-jump logarithm cancellation" begin
  A = ExactComplex[
    1 0 0 1;
    0 -1 2 0;
    0 0 1 -1;
    2 0 0 0
  ]
  B = ExactComplex[
    0 1 0 0;
    1 0 1 0;
    0 -1 2 1;
    0 0 1 -1
  ]
  I4 = Matrix{ExactComplex}(I, 4, 4)
  scalar_no_jump = -2 * I4

  first_channel = A + scalar_no_jump
  second_channel =
    B + scalar_no_jump * A + (1 // 2) * scalar_no_jump * scalar_no_jump
  connected_second = second_channel - (1 // 2) * first_channel * first_channel

  @test connected_second == B - (1 // 2) * A * A
end

@testset "Stinespring TP hierarchy and finite amplitude triangle" begin
  @test static_unitary_jump_tp_coefficient(0) == 1
  for order in 1:8
    @test iszero(static_unitary_jump_tp_coefficient(order))
  end

  for generator_order in 0:6
    for channel_order in 0:(generator_order + 1)
      @test openleg_triangle_contains_channel_order(generator_order, channel_order)
    end
    @test !openleg_triangle_contains_channel_order(
      generator_order, generator_order + 2
    )
  end
end

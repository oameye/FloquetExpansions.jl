using Test
using LinearAlgebra: I, kron

include(joinpath(@__DIR__, "helpers", "physical_gram_channel_reference.jl"))

const RotatingGramExact = Complex{Rational{Int}}
const rotating_gram_im = RotatingGramExact(0 // 1, 1 // 1)

function rotating_gram_recycling_super(
  left::Matrix{RotatingGramExact}, right::Matrix{RotatingGramExact}
)
  return kron(conj.(right), left)
end

function rotating_gram_hamiltonian_super(H::Matrix{RotatingGramExact})
  identity_component = Matrix{RotatingGramExact}(I, 2, 2)
  return -rotating_gram_im *
         (kron(identity_component, H) - kron(transpose(H), identity_component))
end

function rotating_jump_fixture()
  sigma_x = RotatingGramExact[0 1; 1 0]
  sigma_y = RotatingGramExact[0 -rotating_gram_im; rotating_gram_im 0]
  sigma_z = RotatingGramExact[1 0; 0 -1]
  sigma_minus = (1 // 2) * (sigma_x - rotating_gram_im * sigma_y)
  sigma_plus = (1 // 2) * (sigma_x + rotating_gram_im * sigma_y)
  identity_component = Matrix{RotatingGramExact}(I, 2, 2)

  jumps = [
    gram_qr_jump(:loss, 0, sigma_z),
    gram_qr_jump(:loss, 1, rotating_gram_im * sigma_minus),
    gram_qr_jump(:loss, -1, -rotating_gram_im * sigma_plus),
  ]
  # L(t)^dagger L(t) = 2I, so Q = -I after factoring the fixed rate gamma.
  drifts = [gram_qr_drift(0, -identity_component)]
  return (; jumps, drifts, identity_component, sigma_z)
end

@testset "physical Gram orientation is fixed by the Stinespring bra-ket convention" begin
  fixture = rotating_jump_fixture()
  zero_component = zero(fixture.identity_component)
  zero_superoperator = zeros(RotatingGramExact, 4, 4)
  words = gram_stinespring_triangle_words(
    fixture.jumps, fixture.drifts, 1, fixture.identity_component
  )

  @test length(words) == 21

  metric = physical_gram_metric_series(words, 1, zero_component)
  for coefficient in metric, matrix in coefficient.coefficients
    @test matrix == adjoint(matrix)
  end

  channel = physical_gram_channel_series(
    words, 1, zero_superoperator; paste=rotating_gram_recycling_super
  )
  @test formal_period_coefficient(channel[1], 0, zero_superoperator) ==
    Matrix{RotatingGramExact}(I, 4, 4)
end

@testset "physical period-channel logarithm separates Floquet gauge from Van Vleck" begin
  fixture = rotating_jump_fixture()
  zero_superoperator = zeros(RotatingGramExact, 4, 4)
  words = gram_stinespring_triangle_words(
    fixture.jumps, fixture.drifts, 1, fixture.identity_component
  )

  first_order = physical_gram_channel_order_by_phase(
    words, 1, zero_superoperator; paste=rotating_gram_recycling_super
  )
  second_order = physical_gram_channel_order_by_phase(
    words, 2, zero_superoperator; paste=rotating_gram_recycling_super
  )

  @test sort!(collect(keys(first_order))) == [0]
  @test sort!(collect(keys(second_order))) == [-2, -1, 0, 1, 2]

  log_second_order_zero_phase = physical_log_second_order_phase_average(
    first_order[0], second_order[0], zero_superoperator
  )

  @test iszero(
    formal_period_coefficient(log_second_order_zero_phase, 2, zero_superoperator)
  )
  @test iszero(formal_period_coefficient(second_order[-2], 1, zero_superoperator))
  @test iszero(formal_period_coefficient(second_order[2], 1, zero_superoperator))
  @test !iszero(formal_period_coefficient(second_order[-1], 1, zero_superoperator))
  @test !iszero(formal_period_coefficient(second_order[1], 1, zero_superoperator))

  expected = rotating_gram_hamiltonian_super(-(5 // 4) * fixture.sigma_z)
  @test formal_period_coefficient(log_second_order_zero_phase, 1, zero_superoperator) ==
    expected
end

@testset "physical two-output Gram sector resolves the -1 and -1/4 RR pieces" begin
  fixture = rotating_jump_fixture()
  zero_superoperator = zeros(RotatingGramExact, 4, 4)
  words = gram_stinespring_triangle_words(
    fixture.jumps, fixture.drifts, 1, fixture.identity_component
  )
  rr = physical_rr_by_first_mismatch(
    words, zero_superoperator; paste=rotating_gram_recycling_super
  )

  @test sort!(collect(keys(rr))) == [0, 1, 2]
  @test iszero(formal_period_coefficient(rr[0], 1, zero_superoperator))

  m1 = formal_period_coefficient(rr[1], 1, zero_superoperator)
  m2 = formal_period_coefficient(rr[2], 1, zero_superoperator)
  expected_m1 = rotating_gram_hamiltonian_super(-fixture.sigma_z)
  expected_m2 = rotating_gram_hamiltonian_super(-(1 // 4) * fixture.sigma_z)

  @test m1 == expected_m1
  @test m2 == expected_m2
  @test m1 + m2 == rotating_gram_hamiltonian_super(-(5 // 4) * fixture.sigma_z)
end

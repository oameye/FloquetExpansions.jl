using Test
using FloquetExpansions

const CKOutputOverlapExact = Complex{Rational{Int}}
const ck_output_overlap_im = CKOutputOverlapExact(0 // 1, 1 // 1)

function ck_output_model_key(cumulative_harmonics::Vector{Int})
  output_number = length(cumulative_harmonics)
  constraints = FloquetExpansions.CKOutputConstraint{Int}[
    FloquetExpansions.CKOutputConstraint(1, output_stop, harmonic) for
    (output_stop, harmonic) in enumerate(cumulative_harmonics)
  ]
  return FloquetExpansions.CKOutputKernelKey(
    fill(1, output_number),
    [FloquetExpansions.CKOutputBlock(0, output_number, 0, output_number)],
    constraints,
    FloquetExpansions.CKOutputConstraint{Int}[],
  )
end

function ck_output_resolvent_key(poles::Vector{Int})
  return FloquetExpansions.CKOutputKernelKey(
    [1],
    [FloquetExpansions.CKOutputBlock(0, 1, 0, 1)],
    FloquetExpansions.CKOutputConstraint{Int}[],
    FloquetExpansions.CKOutputConstraint{Int}[
      FloquetExpansions.CKOutputConstraint(1, 1, pole) for pole in poles
    ],
  )
end

@testset "finite output overlap reproduces exact one-output Fourier orthogonality" begin
  harmonic_two = ck_output_model_key([2])
  same = FloquetExpansions.ck_output_time_overlap(
    harmonic_two, harmonic_two, ck_output_overlap_im
  )
  @test same.terms == Dict(1 => one(CKOutputOverlapExact))

  harmonic_five = ck_output_model_key([5])
  orthogonal = FloquetExpansions.ck_output_time_overlap(
    harmonic_two, harmonic_five, ck_output_overlap_im
  )
  @test isempty(orthogonal.terms)
end

@testset "finite output overlap reproduces charged physical simplex primitives" begin
  right_two = ck_output_model_key([0, 0])
  left_two = ck_output_model_key([1, 0])
  forward_two = FloquetExpansions.ck_output_time_overlap(
    left_two, right_two, ck_output_overlap_im
  )
  reverse_two = FloquetExpansions.ck_output_time_overlap(
    right_two, left_two, ck_output_overlap_im
  )
  @test forward_two.terms == Dict(1 => ck_output_overlap_im)
  @test reverse_two.terms == Dict(1 => -ck_output_overlap_im)

  right_three = ck_output_model_key([0, 0, 0])
  left_three = ck_output_model_key([2, 1, 0])
  forward_three = FloquetExpansions.ck_output_time_overlap(
    left_three, right_three, ck_output_overlap_im
  )
  @test forward_three.terms == Dict(1 => CKOutputOverlapExact(-1 // 2))

  volume_two = FloquetExpansions.ck_output_time_overlap(
    right_two, right_two, ck_output_overlap_im
  )
  @test volume_two.terms == Dict(2 => CKOutputOverlapExact(1 // 2))
end

@testset "homological and phase pieces pair exactly without sideband enumeration" begin
  free = ck_output_resolvent_key([2])
  phase_five = ck_output_model_key([5])
  overlap = FloquetExpansions.ck_output_time_overlap(free, phase_five, ck_output_overlap_im)
  reverse = FloquetExpansions.ck_output_time_overlap(phase_five, free, ck_output_overlap_im)
  @test overlap.terms == Dict(1 => -ck_output_overlap_im / 3)
  @test reverse.terms == Dict(1 => ck_output_overlap_im / 3)

  excluded = FloquetExpansions.ck_output_time_overlap(
    free, ck_output_model_key([2]), ck_output_overlap_im
  )
  @test isempty(excluded.terms)

  two_poles = ck_output_resolvent_key([2, 5])
  phase_zero = ck_output_model_key([0])
  two_pole_overlap = FloquetExpansions.ck_output_time_overlap(
    two_poles, phase_zero, ck_output_overlap_im
  )
  @test two_pole_overlap.terms == Dict(1 => CKOutputOverlapExact(-1 // 10))

  for excluded_harmonic in (2, 5)
    excluded_overlap = FloquetExpansions.ck_output_time_overlap(
      two_poles, ck_output_model_key([excluded_harmonic]), ck_output_overlap_im
    )
    @test isempty(excluded_overlap.terms)
  end
end

@testset "reverse-bra output overlap is coefficientwise Hermitian" begin
  left = ck_output_resolvent_key([2, 5])
  right = ck_output_model_key([-1])
  forward = FloquetExpansions.ck_output_time_overlap(left, right, ck_output_overlap_im)
  reverse = FloquetExpansions.ck_output_time_overlap(right, left, ck_output_overlap_im)
  @test keys(forward.terms) == keys(reverse.terms)
  for period_power in keys(forward.terms)
    @test forward.terms[period_power] == conj(reverse.terms[period_power])
  end
end

@testset "mismatched physical output channels are orthogonal" begin
  left = ck_output_model_key([0])
  right = FloquetExpansions.CKOutputKernelKey(
    [2],
    [FloquetExpansions.CKOutputBlock(0, 1, 0, 1)],
    [FloquetExpansions.CKOutputConstraint(1, 1, 0)],
    FloquetExpansions.CKOutputConstraint{Int}[],
  )
  overlap = FloquetExpansions.ck_output_time_overlap(left, right, ck_output_overlap_im)
  @test isempty(overlap.terms)
end

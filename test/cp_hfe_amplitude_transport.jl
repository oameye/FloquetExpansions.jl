using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const FE = FloquetExpansions
const SQA = SecondQuantizedAlgebra

space = NLevelSpace(:cp_hfe_transport, 3)
σ11 = Transition(space, :σ, 1, 1)
σ22 = Transition(space, :σ, 2, 2)
σ33 = Transition(space, :σ, 3, 3)
σ12 = Transition(space, :σ, 1, 2)
σ23 = Transition(space, :σ, 2, 3)
σ31 = Transition(space, :σ, 3, 1)
σ13 = Transition(space, :σ, 1, 3)
σ21 = Transition(space, :σ, 2, 1)
σ32 = Transition(space, :σ, 3, 2)
@variables w_cp_transport::Real t_cp_transport::Real γ_cp_transport::Real

w = w_cp_transport
t = t_cp_transport
γ = γ_cp_transport

H0 = 2 * σ11 - σ22 + (3 // 2) * σ33 + σ12 + σ21
H1 = σ12 + 2 * σ23 + im * σ31
H2 = 2 * σ13 - σ21 + im * σ32
H =
  H0 +
  H1 * SQA.expim(-w * t) +
  H1' * SQA.expim(w * t) +
  H2 * SQA.expim(-2 * w * t) +
  H2' * SQA.expim(2 * w * t)

C0 = σ12 + 2 * σ23 + im * σ31
C1 = σ31 + σ12 - im * σ23
Cm1 = 2 * σ21 + σ32
C2 = σ13 - σ21 + im * σ32
C =
  C0 +
  C1 * SQA.expim(-w * t) +
  Cm1 * SQA.expim(w * t) +
  C2 * SQA.expim(-2 * w * t)

liouvillian_vanishes_cp(L::Liouvillian) = iszero(SQA.simplify(L))
periodic_vanishes_cp(G::PeriodicGenerator) =
  all(iszero(SQA.simplify(G[harmonic])) for harmonic in keys(G))

function one_dissipator_cp_reference(H_map::PeriodicGenerator, R::PeriodicGenerator)
  harmonics = filter(!=(0), collect(keys(H_map)))
  H0_map = H_map[0]
  R0 = R[0]

  C_H0 = zero(H0_map)
  C_R0 = zero(H0_map)
  C_3h = zero(H0_map)

  for m in harmonics
    C_H0 -=
      (1 // m^2) * SQA.commutator(SQA.commutator(H_map[m], H0_map), R[-m])
    C_R0 +=
      (1 // (2 * m^2)) * SQA.commutator(H_map[m], SQA.commutator(H_map[-m], R0))
  end

  for m in harmonics, n in harmonics
    if n != m
      C_3h -=
        (1 // (2 * m * n)) *
        SQA.commutator(SQA.commutator(H_map[n], H_map[m - n]), R[-m])
    end
    if m + n != 0
      C_3h -=
        (1 // (2 * m * n)) *
        SQA.commutator(H_map[m], SQA.commutator(H_map[n], R[-(m + n)]))
    end
  end

  return SQA.simplify(C_H0 + C_R0 + C_3h)
end

@testset "native CP-HFE transports physical amplitudes" begin
  reconstruction = FE.cp_hfe_reconstruction(H, w, t, 3, (jump(C, γ),))

  @test length(reconstruction.amplitudes) == 1
  series = only(reconstruction.amplitudes)
  @test series.seed.reference.kind == FE.JUMP_SEED
  @test series.seed.reference.index == 1
  @test series.seed.rate == γ
  @test length(series.coefficients) == 3

  kicks = getfield(reconstruction.coherent, :kick_components)
  A0 = series.seed.amplitude
  expected_A1 = im * SQA.commutator(kicks[1], A0)
  expected_A2 =
    im * SQA.commutator(kicks[2], A0) -
    (1 // 2) * SQA.commutator(kicks[1], SQA.commutator(kicks[1], A0))

  @test periodic_vanishes_cp(series.coefficients[1] - A0)
  @test periodic_vanishes_cp(series.coefficients[2] - expected_A1)
  @test periodic_vanishes_cp(series.coefficients[3] - expected_A2)
end

@testset "finite transported amplitudes reconstruct the averaged Kraus rows" begin
  reconstruction = FE.cp_hfe_reconstruction(H, w, t, 3, (jump(C, γ),))
  series = only(reconstruction.amplitudes)
  finite = FE.finite_transported_amplitude(series)

  direct_average = time_average(harmonics(γ * dissipator(finite(t)), w, t))
  row_average = zero(Liouvillian)
  for channel in reconstruction.channels
    row_average += FE.cp_amplitude_channel_liouvillian(channel)
  end

  @test liouvillian_vanishes_cp(direct_average - row_average)

  expected_generator =
    hamiltonian_action(effective_generator(reconstruction.coherent)) + direct_average
  @test liouvillian_vanishes_cp(reconstruction.generator - expected_generator)
end

@testset "native amplitude transport reproduces the certified one-dissipator sector" begin
  reconstruction = FE.cp_hfe_reconstruction(H, w, t, 3, (jump(C, γ),))
  H_periodic = harmonics(H, w, t)
  H_map = PeriodicGenerator(
    Dict(harmonic => hamiltonian_action(H_periodic[harmonic]) for harmonic in keys(H_periodic)),
    w,
  )
  R = harmonics(γ * dissipator(C), w, t)

  @test liouvillian_vanishes_cp(FE.cp_dissipative_component(reconstruction.amplitudes, 0) - R[0])

  expected_second_order = one_dissipator_cp_reference(H_map, R)
  actual_second_order = FE.cp_dissipative_component(reconstruction.amplitudes, 2)
  @test liouvillian_vanishes_cp(actual_second_order - expected_second_order)
end

@testset "native CP-HFE keeps physical channel provenance" begin
  seeds = FE.physical_amplitude_seeds((collapse(C0), jump(C1, γ)), w, t)

  @test length(seeds) == 2
  @test seeds[1].reference.kind == FE.COLLAPSE_SEED
  @test seeds[1].reference.index == 1
  @test seeds[2].reference.kind == FE.JUMP_SEED
  @test seeds[2].reference.index == 2
end

@testset "periodic scalar jump rates remain outside the first native tranche" begin
  periodic_rate = γ * (1 + cos(w * t))
  channel = jump(C0, periodic_rate)
  @test_throws ArgumentError FE.cp_hfe_reconstruction(H0, w, t, 2, (channel,))
end

using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const FE = FloquetExpansions
const SQA = SecondQuantizedAlgebra

space = NLevelSpace(:native_gksl_transport, 3)
σ11 = Transition(space, :σ, 1, 1)
σ22 = Transition(space, :σ, 2, 2)
σ33 = Transition(space, :σ, 3, 3)
σ12 = Transition(space, :σ, 1, 2)
σ23 = Transition(space, :σ, 2, 3)
σ31 = Transition(space, :σ, 3, 1)
σ13 = Transition(space, :σ, 1, 3)
σ21 = Transition(space, :σ, 2, 1)
σ32 = Transition(space, :σ, 3, 2)
@variables w_native_transport::Real t_native_transport::Real γ_native_transport::Real

w = w_native_transport
t = t_native_transport
γ = γ_native_transport

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
C = C0 + C1 * SQA.expim(-w * t) + Cm1 * SQA.expim(w * t) + C2 * SQA.expim(-2 * w * t)

liouvillian_vanishes(L::Liouvillian) = iszero(FE.canonical_liouvillian(L))
function periodic_vanishes(G::PeriodicGenerator)
  return all(iszero(SQA.simplify(G[harmonic])) for harmonic in keys(G))
end

function one_dissipator_native_reference(H_map::PeriodicGenerator, R::PeriodicGenerator)
  h_harmonics = filter(!=(0), collect(keys(H_map)))
  H0_map = H_map[0]
  R0 = R[0]

  C_H0 = zero(H0_map)
  C_R0 = zero(H0_map)
  C_3h = zero(H0_map)

  for m in h_harmonics
    C_H0 -= (1 // m^2) * SQA.commutator(SQA.commutator(H_map[m], H0_map), R[-m])
    C_R0 += (1 // (2 * m^2)) * SQA.commutator(H_map[m], SQA.commutator(H_map[-m], R0))
  end

  for n in h_harmonics, k in h_harmonics
    generated = n + k
    iszero(generated) && continue
    C_3h -=
      (1 // (2 * generated * n)) *
      SQA.commutator(SQA.commutator(H_map[n], H_map[k]), R[-generated])
  end

  for m in h_harmonics, n in h_harmonics
    if m + n != 0
      C_3h -=
        (1 // (2 * m * n)) * SQA.commutator(H_map[m], SQA.commutator(H_map[n], R[-(m + n)]))
    end
  end

  return SQA.simplify(C_H0 + C_R0 + C_3h)
end

@testset "native GKSL expansion transports jump amplitudes by the coherent micromotion" begin
  native = FE.native_gksl_expansion(H, w, t, 3, (jump(C, γ),))

  @test length(native.amplitudes) == 1
  transported = only(native.amplitudes)
  @test transported.seed.reference.kind == FE.JUMP_SEED
  @test transported.seed.reference.index == 1
  @test transported.seed.rate == γ
  @test length(transported.coefficients) == 3

  K = native.coherent.micromotion_components
  A0 = transported.seed.amplitude
  expected_A1 = im * SQA.commutator(K[1], A0)
  expected_A2 =
    im * SQA.commutator(K[2], A0) -
    (1 // 2) * SQA.commutator(K[1], SQA.commutator(K[1], A0))

  @test periodic_vanishes(transported.coefficients[1] - A0)
  @test periodic_vanishes(transported.coefficients[2] - expected_A1)
  @test periodic_vanishes(transported.coefficients[3] - expected_A2)
end

@testset "harmonic jump channels reproduce the period-averaged dissipator" begin
  native = FE.native_gksl_expansion(H, w, t, 3, (jump(C, γ),))
  finite = FE.finite_transported_amplitude(only(native.amplitudes))

  direct_average = time_average(harmonics(γ * dissipator(finite(t)), w, t))
  channel_sum = zero(Liouvillian)
  for channel in native.channels
    channel_sum += FE.channel_liouvillian(channel)
  end

  @test liouvillian_vanishes(direct_average - channel_sum)

  expected_generator =
    hamiltonian_action(effective_generator(native.coherent)) + direct_average
  @test liouvillian_vanishes(native.generator - expected_generator)
end

@testset "native GKSL expansion reproduces the certified one-dissipator sector" begin
  native = FE.native_gksl_expansion(H, w, t, 3, (jump(C, γ),))
  H_periodic = harmonics(H, w, t)
  H_map = PeriodicGenerator(
    Dict(
      harmonic => hamiltonian_action(H_periodic[harmonic]) for harmonic in keys(H_periodic)
    ),
    w,
  )
  R = harmonics(γ * dissipator(C), w, t)

  @test !liouvillian_vanishes(R[3])
  @test !liouvillian_vanishes(R[-3])
  @test liouvillian_vanishes(FE.native_dissipative_component(native.amplitudes, 0) - R[0])

  expected_second_order = one_dissipator_native_reference(H_map, R)
  actual_second_order = FE.native_dissipative_component(native.amplitudes, 2)
  @test liouvillian_vanishes(actual_second_order - expected_second_order)
end

@testset "native GKSL expansion keeps microscopic channel provenance" begin
  seeds = FE.jump_amplitude_seeds((collapse(C0), jump(C1, γ), collapse(C2)), w, t)

  @test [seed.reference.kind for seed in seeds] == [FE.COLLAPSE_SEED, FE.JUMP_SEED, FE.COLLAPSE_SEED]
  @test [seed.reference.index for seed in seeds] == [1, 1, 2]
  @test seeds[3].amplitude == harmonics(C2, w, t)
end

@testset "native GKSL expansion rejects periodic scalar jump rates" begin
  channel = jump(C0, γ * (1 + cos(w * t)))
  failure = @test_throws ArgumentError FE.native_gksl_expansion(H0, w, t, 2, (channel,))
  @test occursin("time-independent jump rate", sprint(showerror, failure.value))
end

using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const FE = FloquetExpansions
const SQA = SecondQuantizedAlgebra

liouvillian_vanishes(L::Liouvillian) = iszero(FE.canonical_liouvillian(L))

function static_similarity_generator(H_map::PeriodicGenerator, R::PeriodicGenerator)
  result = zero(R[0])
  for harmonic in keys(H_map)
    iszero(harmonic) && continue
    result += (1 // (2 * harmonic^2)) * SQA.commutator(H_map[harmonic], R[-harmonic])
  end
  return result
end

# Polarize the direct Liouvillian van Vleck expansion of H + λR in a bookkeeping scale λ on the
# dissipator. Its order-1 effective component is quadratic and its order-2 component cubic in λ,
# so these stencils return exact polynomial coefficients: the two-dissipator sector at first order
# and the one-dissipator sector at second order.
function direct_dissipator_sectors(H_map::PeriodicGenerator, R::PeriodicGenerator, ω)
  first_order = Dict{Int,Liouvillian}()
  second_order = Dict{Int,Liouvillian}()
  for λ in (-1, 0, 1, 2)
    expansion = floquet_expansion(H_map + λ * R, VanVleck(), 3)
    first_order[λ] = ω * effective_component(expansion, 1)
    second_order[λ] = ω^2 * effective_component(expansion, 2)
  end

  two_dissipator_first = (1 // 2) * (first_order[1] + first_order[-1]) - first_order[0]
  one_dissipator_second =
    (1 // 6) *
    (6 * second_order[1] - 2 * second_order[-1] - 3 * second_order[0] - second_order[2])
  return (; first=first_order[1], two_dissipator_first, one_dissipator_second)
end

function validate_periodic_jump_sectors(H, channel, ω, t)
  native = FE.native_gksl_expansion(H, ω, t, 3, (channel,))
  H_map = harmonics(hamiltonian_action(H), ω, t)
  R = harmonics(liouvillian(zero(SQA.QAdd); channels=(channel,)), ω, t)
  direct = direct_dissipator_sectors(H_map, R, ω)

  # At first order the native construction reproduces the coherent and one-dissipator sectors of
  # the direct expansion exactly. The remainder is the two-dissipator sector of order γ²/ω, which
  # lies outside the native GKSL expansion and is nonzero for a periodic jump.
  retained_first =
    hamiltonian_action(ω * effective_component(native.coherent, 1)) +
    FE.native_dissipative_component(native.amplitudes, 1)
  @test !liouvillian_vanishes(direct.two_dissipator_first)
  @test liouvillian_vanishes(direct.first - retained_first - direct.two_dissipator_first)

  # At second order the one-dissipator sector of the direct expansion and the native coefficient
  # are independently obtained representatives that differ exactly by the static similarity
  # [B_R^(2), H_0], which is nonzero here.
  similarity = SQA.commutator(static_similarity_generator(H_map, R), H_map[0])
  native_second = FE.native_dissipative_component(native.amplitudes, 2)
  @test !liouvillian_vanishes(similarity)
  @test liouvillian_vanishes(direct.one_dissipator_second - native_second - similarity)
  return nothing
end

@testset "native GKSL expansion: sideband jump separates the one- and two-dissipator sectors" begin
  pauli = PauliSpace(:native_gksl_sideband_jump)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  σminus = (σx - im * σy) / 2

  @variables ω_sideband_jump::Real t_sideband_jump::Real Ω_sideband_jump::Real
  @variables Δ_sideband_jump::Real ε_sideband_jump::Real γ_sideband_jump::Real
  ω = ω_sideband_jump
  t = t_sideband_jump
  Ω = Ω_sideband_jump
  Δ = Δ_sideband_jump
  ε = ε_sideband_jump
  γ = γ_sideband_jump

  # L(t) = σ₋ + ε e^{-iωt} σz carries the harmonics 0 and +1, so its dissipator has the harmonics
  # -1, 0 and 1. Every parameter stays symbolic.
  sideband = cos(ω * t) - im * sin(ω * t)
  L = σminus + ε * sideband * σz
  H = Δ * σz + Ω * cos(ω * t) * σx

  validate_periodic_jump_sectors(H, jump(L, γ), ω, t)
end

@testset "native GKSL expansion: rotating jump separates the one- and two-dissipator sectors" begin
  pauli = PauliSpace(:native_gksl_rotating_jump)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω_rotating_jump::Real t_rotating_jump::Real
  ω = ω_rotating_jump
  t = t_rotating_jump

  # Exact nondegenerate point Δ = 2, Ω = 1, η = 1 of
  # H(t) = Δ σz / 2 + Ω cos(ωt) σx and
  # J(t) = σz + η [cos(ωt) σy + sin(ωt) σx].
  H = σz + cos(ω * t) * σx
  rotating_jump = σz + cos(ω * t) * σy + sin(ω * t) * σx

  # Directional modulation, rather than a scalar rate modulation, produces genuinely distinct
  # dissipative harmonics. The dissipator of this jump has the harmonics -2 to 2, so the sectors
  # also couple through the second sidebands, which the sideband jump lacks.
  validate_periodic_jump_sectors(H, jump(rotating_jump, 1), ω, t)
end

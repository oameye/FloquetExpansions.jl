using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const FE = FloquetExpansions
const SQA = SecondQuantizedAlgebra

periodic_jump_zero(L::Liouvillian) = iszero(FE.canonical_liouvillian(L))

function periodic_jump_retained_component(reconstruction, grade::Int)
  coherent = getfield(reconstruction.coherent, :effective_components)[grade + 1]
  dissipative = FE.cp_dissipative_component(reconstruction.amplitudes, grade)
  return SQA.simplify(hamiltonian_action(coherent) + dissipative)::Liouvillian
end

function periodic_jump_substitute(L::Liouvillian, rules)
  result = zero(L)
  for (left, right, coefficient) in terms(L)
    result += FE.action(left, right, SQA.substitute_cnum(coefficient, rules))
  end
  return FE.canonical_liouvillian(result)
end

function periodic_jump_linear_rate(L::Liouvillian, rate)
  values = ntuple(
    n -> periodic_jump_substitute(L, Dict(rate => n - 1)),
    4,
  )
  Δ1 = values[2] - values[1]
  Δ2 = values[3] - 2 * values[2] + values[1]
  Δ3 = values[4] - 3 * values[3] + 3 * values[2] - values[1]
  return FE.canonical_liouvillian(Δ1 - (1 // 2) * Δ2 + (1 // 3) * Δ3)
end

@testset "physical periodic jump separates one-R gauge physics from RR" begin
  pauli = PauliSpace(:cp_hfe_periodic_jump)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  σminus = (σx - im * σy) / 2

  @variables ω_cp_jump::Real t_cp_jump::Real Ω_cp_jump::Real Δ_cp_jump::Real
  @variables ε_cp_jump::Real γ_cp_jump::Real
  ω = ω_cp_jump
  t = t_cp_jump
  Ω = Ω_cp_jump
  Δ = Δ_cp_jump
  ε = ε_cp_jump
  γ = γ_cp_jump

  sideband = cos(ω * t) - im * sin(ω * t)
  L = σminus + ε * sideband * σz
  H = Δ * σz + Ω * cos(ω * t) * σx
  channels = (jump(L, γ),)

  R = harmonics(liouvillian(0 * σx; channels=channels), ω, t)
  H_harmonics = harmonics(H, ω, t)
  H_map = PeriodicGenerator(
    Dict(harmonic => hamiltonian_action(H_harmonics[harmonic]) for harmonic in keys(H_harmonics)),
    ω,
  )

  h_harmonics = filter(!=(0), collect(keys(H_map)))
  r_harmonics = filter(!=(0), collect(keys(R)))
  H0 = H_map[0]
  R0 = R[0]

  # The microscopic periodic jump has a genuine one-dissipator static frame correction.
  B_R = zero(R0)
  for m in h_harmonics
    B_R += (1 // (2 * m^2)) * SQA.commutator(H_map[m], R[-m])
  end
  @test !periodic_jump_zero(B_R)

  # At first inverse-frequency order the direct Liouvillian expansion contains a genuine RR
  # term. Native amplitude transport reproduces the coherent + one-R content, so the exact
  # difference is the quadratic dissipative sector.
  R_osc = PeriodicGenerator(Dict(m => R[m] for m in r_harmonics), ω)
  RR = zero(R0)
  for m in r_harmonics
    RR += (im // m) * compose(R_osc[-m], R_osc[m])
  end
  @test !periodic_jump_zero(RR)

  native = FE.cp_hfe_reconstruction(H, ω, t, 3, channels)
  raw = floquet_expansion(H, ω, t, VanVleck(), 3; channels=channels)
  raw_components = getfield(raw, :effective_components)
  @test periodic_jump_zero(
    raw_components[2] - periodic_jump_retained_component(native, 1) - RR
  )

  # A second-order Van Vleck coefficient contains at most three dissipative vertices, hence is a
  # cubic polynomial in the formal channel strength γ. Exact forward differences at γ=0,1,2,3
  # isolate its coefficient linear in γ without reusing the #104 one-R formula. The native
  # dissipative coefficient is itself linear in γ. Their independently obtained representatives
  # differ exactly by the certified static similarity [B_R^(2), H_0].
  direct_one_R2 = periodic_jump_linear_rate(raw_components[3], γ)
  native_R2 = FE.cp_dissipative_component(native.amplitudes, 2)
  native_one_R2 = periodic_jump_substitute(native_R2, Dict(γ => 1))
  gauge_one_R2 = periodic_jump_substitute(SQA.commutator(B_R, H0), Dict(γ => 1))
  @test periodic_jump_zero(direct_one_R2 - native_one_R2 - gauge_one_R2)

  # With the Hamiltonian switched off there is no coherent amplitude transport at O(1/ω), while
  # the same microscopic periodic jump still has the nonzero direct-Liouvillian RR contribution.
  dissipative_only = FE.cp_hfe_reconstruction(0 * σx, ω, t, 2, channels)
  @test iszero(only(dissipative_only.amplitudes).coefficients[2])
end

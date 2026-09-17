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

  # Build the complete one-R second-order sector independently from the #104 identities. This
  # isolates the theorem from the RR/RRR pieces carried by the full direct Liouvillian expansion.
  C_H0 = zero(R0)
  C_R0 = zero(R0)
  V_R_a = zero(R0)
  for m in h_harmonics
    C_H0 -= (1 // m^2) * SQA.commutator(SQA.commutator(H_map[m], H0), R[-m])
    C_R0 +=
      (1 // (2 * m^2)) * SQA.commutator(H_map[m], SQA.commutator(H_map[-m], R0))
    V_R_a -=
      (1 // (2 * m^2)) * (
        SQA.commutator(R[-m], SQA.commutator(H0, H_map[m])) +
        SQA.commutator(H_map[-m], SQA.commutator(R0, H_map[m])) +
        SQA.commutator(H_map[-m], SQA.commutator(H0, R[m]))
      )
  end

  three_harmonic_bound = max(2 * maximum(abs, h_harmonics), maximum(abs, r_harmonics))
  three_harmonics = filter(!=(0), collect((-three_harmonic_bound):three_harmonic_bound))
  C_3h = zero(R0)
  V_R_b = zero(R0)
  for m in three_harmonics, n in three_harmonics
    if n != m
      C_3h -=
        (1 // (2 * m * n)) *
        SQA.commutator(SQA.commutator(H_map[n], H_map[m - n]), R[-m])
      V_R_b -=
        (1 // (3 * n * m)) * (
          SQA.commutator(R[-n], SQA.commutator(H_map[n - m], H_map[m])) +
          SQA.commutator(H_map[-n], SQA.commutator(R[n - m], H_map[m])) +
          SQA.commutator(H_map[-n], SQA.commutator(H_map[n - m], R[m]))
        )
    end
    if m + n != 0
      C_3h -=
        (1 // (2 * m * n)) *
        SQA.commutator(H_map[m], SQA.commutator(H_map[n], R[-(m + n)]))
    end
  end

  R_CP = SQA.simplify(C_H0 + C_R0 + C_3h)
  V_R = SQA.simplify(V_R_a + V_R_b)
  native_R2 = FE.cp_dissipative_component(native.amplitudes, 2)

  @test periodic_jump_zero(native_R2 - R_CP)
  @test periodic_jump_zero(V_R - native_R2 - SQA.commutator(B_R, H0))

  # With the Hamiltonian switched off there is no coherent amplitude transport at O(1/ω), while
  # the same microscopic periodic jump still has the nonzero direct-Liouvillian RR contribution.
  dissipative_only = FE.cp_hfe_reconstruction(0 * σx, ω, t, 2, channels)
  @test iszero(only(dissipative_only.amplitudes).coefficients[2])
end

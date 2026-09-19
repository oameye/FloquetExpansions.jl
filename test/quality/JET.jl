using Test
using FloquetExpansions
using JET: JET
using Symbolics: @variables

@testset "JET" begin
  JET.test_package(FloquetExpansions; target_modules=(FloquetExpansions,))
end

@testset "graded output optimizer stability" begin
  left = FloquetExpansions.graded_output_series(Int, 0 // 1, 3)
  right = FloquetExpansions.graded_output_series(Int, 0 // 1, 3)
  FloquetExpansions.graded_output_accumulate!(left, 1, 1, 2 // 1)
  FloquetExpansions.graded_output_accumulate!(right, 1, 2, 3 // 1)
  compose_output =
    (left_key, right_key) -> FloquetExpansions.OutputSectorComposition(left_key + right_key)
  algebra = FloquetExpansions.GradedOutputProduct(compose_output, *)

  JET.@test_opt target_modules=(FloquetExpansions,) algebra(left, right)
end

@testset "completion optimizer stability" begin
  fock = FockSpace(:jet_completion_fock)
  a = Destroy(fock, :a)
  @variables ω::Real t::Real
  gram_frame = DissipativeFrame(a, a^2)
  gram_generator = liouvillian(0 * a; channels=(collapse(a + a^2), collapse(a + im * a^2)))
  gram_expansion = floquet_expansion(gram_generator, ω, t, VanVleck(), 1)

  JET.@test_opt target_modules=(FloquetExpansions,) positive_completion(
    gram_expansion, Gram(), gram_frame
  )

  pauli = PauliSpace(:jet_completion_pauli)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  recursive_frame = DissipativeFrame(σx, σy, σz)
  recursive_expansion = floquet_expansion(
    cos(ω * t) * σx, ω, t, VanVleck(), 3; channels=(collapse(σz),)
  )

  JET.@test_opt target_modules=(FloquetExpansions,) positive_completion(
    recursive_expansion, Gram(), recursive_frame
  )

  spectral_frame = DissipativeFrame(σz, σy)
  JET.@test_opt target_modules=(FloquetExpansions,) positive_completion(
    recursive_expansion, Spectral(), spectral_frame
  )

  completion = positive_completion(gram_expansion, Gram(), gram_frame)
  JET.@test_opt target_modules=(FloquetExpansions,) hamiltonian(completion)
  JET.@test_opt target_modules=(FloquetExpansions,) kossakowski_component(completion, 0)
end

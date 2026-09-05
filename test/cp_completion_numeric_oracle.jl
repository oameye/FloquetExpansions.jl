using Test
using FloquetExpansions
using LinearAlgebra: eigvals
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables

const SQA = SecondQuantizedAlgebra

function numeric_kossakowski(matrix, substitutions)
  result = zeros(ComplexF64, size(matrix))
  for index in eachindex(matrix)
    value = SQA.to_num(matrix[index])
    re = Symbolics.value(Symbolics.substitute(real(value), substitutions))
    im = Symbolics.value(Symbolics.substitute(imag(value), substitutions))
    result[index] = complex(Float64(re), Float64(im))
  end
  return result
end

function minimum_kossakowski_eigenvalue(matrix, substitutions)
  values = eigvals(numeric_kossakowski(matrix, substitutions))
  @test maximum(abs.(imag.(values))) <= 1.0e-10
  return minimum(real.(values))
end

@testset "completion removes the weak negative eigenvalue of a rank-one truncation" begin
  pauli = PauliSpace(:cp_numeric_weak_channel)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  @variables ω::Real

  coherent_rotation = hamiltonian_action(σx)
  diagonal_difference = dissipator(σy) - dissipator(σz)
  generator = PeriodicGenerator(
    Dict(0 => dissipator(σz), 1 => coherent_rotation, -1 => diagonal_difference), ω
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)
  gram = positive_completion(expansion, Gram(), frame)
  spectral = positive_completion(expansion, Spectral(), frame)
  substitutions = Dict(ω => 10.0)

  raw_minimum = minimum_kossakowski_eigenvalue(kossakowski(expansion, frame), substitutions)
  gram_minimum = minimum_kossakowski_eigenvalue(kossakowski(gram), substitutions)
  spectral_minimum = minimum_kossakowski_eigenvalue(kossakowski(spectral), substitutions)

  @test raw_minimum < -1.0e-8
  @test gram_minimum >= -1.0e-10
  @test spectral_minimum >= -1.0e-10
end

@testset "driven-qubit completions are numerically positive in the adapted frame" begin
  pauli = PauliSpace(:cp_numeric_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  bright = σy + σz
  dark = σy - σz
  frame = DissipativeFrame(bright, dark)
  @variables ω::Real t::Real Ω::Real

  expansion = floquet_expansion(
    Ω * cos(ω * t) * σx, ω, t, VanVleck(), 3; channels=(collapse(bright),)
  )
  gram = positive_completion(expansion, Gram(), frame)
  spectral = positive_completion(expansion, Spectral(), frame)
  substitutions = Dict(Ω => 0.7, ω => 12.0)

  @test minimum_kossakowski_eigenvalue(kossakowski(gram), substitutions) >= -1.0e-10
  @test minimum_kossakowski_eigenvalue(kossakowski(spectral), substitutions) >= -1.0e-10
end

@testset "bosonic modulated-loss completion is positive without a Fock cutoff" begin
  fock = FockSpace(:cp_numeric_modulated_loss)
  a = Destroy(fock, :a)
  frame = DissipativeFrame(a, a^2, a' * a^2)
  @variables ω::Real t::Real K::Real

  H = K * cos(ω * t) * a'^2 * a^2
  expansion = floquet_expansion(
    H,
    ω,
    t,
    VanVleck(),
    2;
    channels=(jump(a, 2 + sin(ω * t)), collapse(a^2)),
  )
  completion = positive_completion(expansion, Gram(), frame)
  substitutions = Dict(K => 0.6, ω => 12.0)

  raw_minimum = minimum_kossakowski_eigenvalue(kossakowski(expansion, frame), substitutions)
  completed_minimum = minimum_kossakowski_eigenvalue(kossakowski(completion), substitutions)

  @test raw_minimum < -1.0e-10
  @test completed_minimum >= -1.0e-10
end

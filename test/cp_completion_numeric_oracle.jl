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

@testset "completion removes a weak negative eigenvalue of a rank-one truncation" begin
  pauli = PauliSpace(:cp_numeric_weak_channel)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  frame = DissipativeFrame(σz, σy)
  @variables ω::Real

  coherent_rotation = hamiltonian_action(σx)
  diagonal_difference = dissipator(σy) - dissipator(σz)
  quadrature = im * (1 // 2) * diagonal_difference
  generator = PeriodicGenerator(
    Dict(
      0 => dissipator(σz),
      1 => coherent_rotation + quadrature,
      -1 => coherent_rotation - quadrature,
    ),
    ω,
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

@testset "bosonic Gram closure is positive without a Fock cutoff" begin
  fock = FockSpace(:cp_numeric_modulated_loss)
  a = Destroy(fock, :a)
  frame = DissipativeFrame(a, a^2, a' * a^2)
  @variables ω::Real

  D1 = dissipator(a)
  D2 = dissipator(a^2)
  quadrature = im * (1 // 2) * D2
  generator = PeriodicGenerator(
    Dict(0 => D1 + D2, 1 => (1 // 2) * D1 + quadrature, -1 => (1 // 2) * D1 - quadrature), ω
  )
  expansion = floquet_expansion(generator, VanVleck(), 2)
  completion = positive_completion(expansion, Gram(), frame)
  substitutions = Dict(ω => 12.0)

  raw_minimum = minimum_kossakowski_eigenvalue(kossakowski(expansion, frame), substitutions)
  completed_minimum = minimum_kossakowski_eigenvalue(kossakowski(completion), substitutions)

  @test raw_minimum < -1.0e-10
  @test completed_minimum >= -1.0e-10
end

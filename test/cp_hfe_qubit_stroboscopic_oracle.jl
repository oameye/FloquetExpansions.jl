using Test
using FloquetExpansions
using LinearAlgebra: I, kron, tr
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables

const SQA = SecondQuantizedAlgebra

function cp_oracle_scalar(value, substitutions)
  numeric = SQA.to_num(value)
  scalar_rules = Dict(
    Symbolics.unwrap(parameter) => replacement for (parameter, replacement) in substitutions
  )
  re = Symbolics.value(Symbolics.substitute(Symbolics.unwrap(real(numeric)), scalar_rules))
  im = Symbolics.value(Symbolics.substitute(Symbolics.unwrap(imag(numeric)), scalar_rules))
  return complex(Float64(re), Float64(im))
end

function cp_oracle_pauli_matrix(operator::SQA.QAdd, substitutions)
  σx = ComplexF64[0 1; 1 0]
  σy = ComplexF64[0 -im; im 0]
  σz = ComplexF64[1 0; 0 -1]
  paulis = (σx, σy, σz)
  result = zeros(ComplexF64, 2, 2)
  identity = Matrix{ComplexF64}(I, 2, 2)

  for (term, coefficient) in operator
    term_matrix = identity
    for op in term.ops
      term_matrix = term_matrix * paulis[Int(op.l1)]
    end
    result .+= cp_oracle_scalar(coefficient, substitutions) .* term_matrix
  end
  return result
end

function cp_oracle_liouvillian_matrix(generator::Liouvillian, substitutions)
  result = zeros(ComplexF64, 4, 4)
  for (left, right, coefficient) in terms(generator)
    left_matrix = cp_oracle_pauli_matrix(left, substitutions)
    right_matrix = cp_oracle_pauli_matrix(right, substitutions)
    scalar = cp_oracle_scalar(coefficient, substitutions)
    result .+= scalar .* kron(transpose(right_matrix), left_matrix)
  end
  return result
end

function cp_oracle_dissipator_matrix(collapse_operator::Matrix{ComplexF64})
  identity = Matrix{ComplexF64}(I, 2, 2)
  norm = collapse_operator' * collapse_operator
  return kron(conj(collapse_operator), collapse_operator) -
         (1 // 2) * (kron(identity, norm) + kron(transpose(norm), identity))
end

function cp_oracle_exact_qubit_map(ω::Float64, Ω::Float64; steps::Int=2048)
  σx = ComplexF64[0 1; 1 0]
  σy = ComplexF64[0 -im; im 0]
  σz = ComplexF64[1 0; 0 -1]
  identity = Matrix{ComplexF64}(I, 2, 2)
  bright = σy + σz
  dissipative = cp_oracle_dissipator_matrix(bright)

  period = 2π / ω
  dt = period / steps
  propagator = Matrix{ComplexF64}(I, 4, 4)
  for step in 1:steps
    phase = 2π * (step - 1 // 2) / steps
    hamiltonian = Ω * cos(phase) * σx
    coherent = -im * (kron(identity, hamiltonian) - kron(transpose(hamiltonian), identity))
    propagator = exp(dt * (coherent + dissipative)) * propagator
  end
  return propagator
end

function cp_oracle_trace_signature(propagator::Matrix{ComplexF64})
  return [tr(propagator^power) / 4 for power in 1:4]
end

function cp_oracle_trace_error(exact, approximate)
  exact_signature = cp_oracle_trace_signature(exact)
  approximate_signature = cp_oracle_trace_signature(approximate)
  return maximum(abs.(exact_signature - approximate_signature))
end

@testset "driven-qubit Floquet spectrum converges to independent one-period dynamics" begin
  pauli = PauliSpace(:cp_hfe_stroboscopic_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σy = Pauli(pauli, :sigma, 2)
  σz = Pauli(pauli, :sigma, 3)
  bright = σy + σz
  dark = σy - σz
  adapted = DissipativeFrame(bright, dark)
  @variables ω_cp_oracle::Real t_cp_oracle::Real Ω_cp_oracle::Real

  ω = ω_cp_oracle
  t = t_cp_oracle
  Ω = Ω_cp_oracle
  H = Ω * cos(ω * t) * σx
  channels = (collapse(bright),)

  native_order1 = FloquetExpansions.cp_hfe_reconstruction(H, ω, t, 1, channels)
  native_order3 = FloquetExpansions.cp_hfe_reconstruction(H, ω, t, 3, channels)
  raw = floquet_expansion(H, ω, t, VanVleck(), 3; channels=channels)
  gram = positive_completion(raw, Gram(), adapted)
  spectral = positive_completion(raw, Spectral(), adapted)

  generators = (
    native_order1=native_order1.generator,
    native_order3=native_order3.generator,
    raw_order3=effective_generator(raw),
    gram_order3=effective_generator(gram),
    spectral_order3=effective_generator(spectral),
  )

  frequencies = (8.0, 12.0, 18.0)
  drive = 0.7
  errors = Dict(name => Float64[] for name in keys(generators))

  for frequency in frequencies
    exact_map = cp_oracle_exact_qubit_map(frequency, drive)
    period = 2π / frequency
    substitutions = Dict(ω => frequency, Ω => drive)

    for name in keys(generators)
      effective_matrix = cp_oracle_liouvillian_matrix(generators[name], substitutions)
      approximate_map = exp(period * effective_matrix)
      push!(errors[name], cp_oracle_trace_error(exact_map, approximate_map))
    end
  end

  for method_errors in values(errors)
    @test all(isfinite, method_errors)
    @test last(method_errors) < first(method_errors)
  end

  # The second-order physical amplitude transport must improve the high-frequency stroboscopic
  # spectrum over the undressed leading GKSL generator on this same Gram/Spectral benchmark.
  @test last(errors[:native_order3]) < last(errors[:native_order1])
end

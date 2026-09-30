using Test
using FloquetExpansions
using LinearAlgebra: eigvals
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables

space = FockSpace(:resonator)
a = Destroy(space, :a)
@variables ω₀::Real ω₁::Real ω₂::Real ω̄::Real K::Real κ::Real n̄::Real t::Real

function numeric(coefficient, values::AbstractDict)
  z = SQA.to_num(coefficient)
  part(x) = Float64(Symbolics.value(Symbolics.substitute(x, values)))
  return complex(part(real(z)), part(imag(z)))
end

heisenberg(L::Liouvillian, O) = sum(c * B * O * A for (A, B, c) in terms(L))

function first_moment_matrix(evolve, modes, values)
  basis = vcat(modes, adjoint.(modes))
  M = zeros(ComplexF64, length(basis), length(basis))
  for (row, operator) in enumerate(basis)
    for (term, coefficient) in evolve(operator)
      column = findfirst(b -> [b] == term.ops, basis)
      M[row, column] += numeric(coefficient, values)
    end
  end
  return M
end

function stationary_quadratic_moments(evolve, modes, values)
  basis = vcat(modes, adjoint.(modes))
  monomials = unique(
    vcat(
      [
        [term.ops for (term, _) in u * v if length(term.ops) == 2] for u in basis for
        v in basis
      ]...,
    ),
  )
  A = zeros(ComplexF64, length(monomials), length(monomials))
  b = zeros(ComplexF64, length(monomials))
  for (row, ops) in enumerate(monomials)
    for (term, coefficient) in evolve(prod(ops))
      if isempty(term.ops)
        b[row] -= numeric(coefficient, values)
      else
        A[row, findfirst(==(term.ops), monomials)] += numeric(coefficient, values)
      end
    end
  end
  return Dict(zip(monomials, A \ b))
end

function quadratic_expectation(operator, moments, values)
  return sum(
    (isempty(term.ops) ? 1 : moments[term.ops]) * numeric(coefficient, values) for
    (term, coefficient) in substitute(operator, values)
  )
end

@testset "the carrier embedding preserves the physical algebra" begin
  embedding = CarrierEmbedding(a, [ω₁, (ω₁ + ω₂) / 2, ω₂]; t)
  image = reconstruct(embedding, a)
  @test isequal(commutator(image, image'), one(image))
  @test isequal(
    reconstruct(embedding, a * a' * a), image * reconstruct(embedding, a') * image
  )
end

@testset "number-conserving linear dynamics reconstructs exactly" begin
  embedding = CarrierEmbedding(a, [ω₁, ω₂]; t)
  H = ω₀ * a' * a
  carrier_H = harmonic_balance(H, embedding)
  image = reconstruct(embedding, a)
  explicit_derivative = sum(
    Symbolics.derivative(coefficient, t) * prod(term.ops) for
    (term, coefficient) in exponential_form(image)
  )
  defect =
    explicit_derivative + im * commutator(carrier_H, image) -
    reconstruct(embedding, im * commutator(H, a))
  @test iszero(simplify(defect))
end

@testset "damped counter-rotating dynamics keeps the physical poles" begin
  embedding = CarrierEmbedding(a, [ω₁, ω₂]; t)
  x² = (a + a')^2 / 2
  p² = -(a - a')^2 / 2
  H = p² / 2 + ω₀^2 * x² / 2
  L = liouvillian(H; channels=(jump(a, κ * (n̄ + 1)), jump(a', κ * n̄)))
  carrier_L = harmonic_balance(L, embedding)
  values = Dict(ω₀ => 1.3, ω₁ => 0.9, ω₂ => 1.7, κ => 0.2, n̄ => 0.5)

  physical = eigvals(first_moment_matrix(O -> heisenberg(L, O), [a], values))
  modes = collect(Iterators.flatten(carrier_modes(embedding)))
  carrier = eigvals(first_moment_matrix(O -> heisenberg(carrier_L, O), modes, values))
  expected = [z + s * im * values[ω] for z in physical for s in (-1, 1) for ω in (ω₁, ω₂)]
  by = z -> (round(imag(z); digits=8), round(real(z); digits=8))
  @test sort(carrier; by) ≈ sort(expected; by)

  physical_moments = stationary_quadratic_moments(O -> heisenberg(L, O), [a], values)
  carrier_moments = stationary_quadratic_moments(
    O -> heisenberg(carrier_L, O), modes, values
  )
  for observable in (a' * a, a * a)
    exact = quadratic_expectation(observable, physical_moments, values)
    carrier_observable = reconstruct(embedding, observable)
    for time in (0.0, 0.7, 2.3)
      at_time = merge(values, Dict(t => time))
      @test quadratic_expectation(carrier_observable, carrier_moments, at_time) ≈ exact
    end
  end
end

@testset "resonances follow the declared carrier frequencies" begin
  H = K * a' * a' * a * a
  midpoint = CarrierEmbedding(a, [ω₁, (ω₁ + ω₂) / 2, ω₂]; t)
  (one₊, _), (middle₊, _), (two₊, _) = carrier_modes(midpoint)
  conversion = only(collect(one₊' * two₊' * middle₊ * middle₊)).first.ops
  coefficient(H, ops) = [c for (term, c) in H if term.ops == ops]
  @test isequal(coefficient(harmonic_balance(H, midpoint), conversion), [K / 3])

  independent = CarrierEmbedding(a, [ω₁, ω̄, ω₂]; t)
  (one₊, _), (middle₊, _), (two₊, _) = carrier_modes(independent)
  conversion = only(collect(one₊' * two₊' * middle₊ * middle₊)).first.ops
  @test isempty(coefficient(harmonic_balance(H, independent), conversion))
end

@testset "invalid embeddings and generators are rejected" begin
  @test_throws ArgumentError CarrierEmbedding(a, [0.3 * ω₁]; t)
  @test_throws ArgumentError CarrierEmbedding(a, [ω₁ * t]; t)
  @test_throws ArgumentError CarrierEmbedding(a', [ω₁]; t)
  embedding = CarrierEmbedding(a, [ω₁]; t)
  other = Destroy(FockSpace(:other), :b)
  @test_throws ArgumentError reconstruct(embedding, other)
  @test_throws ArgumentError harmonic_balance(t * a' * a, embedding)
end

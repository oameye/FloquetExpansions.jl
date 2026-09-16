using Test
using FloquetExpansions
using LinearAlgebra: norm
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

include(joinpath(@__DIR__, "helpers", "shared.jl"))
include(joinpath(@__DIR__, "helpers", "bloch_feshbach_reference.jl"))

liouvillian_vanishes(L::Liouvillian) = iszero(SQA.simplify(L))
function liouvillian_vanishes(G::PeriodicGenerator{Liouvillian})
  return all(liouvillian_vanishes(G[harmonic]) for harmonic in keys(G))
end

hamiltonian_product(left::SQA.QAdd, right::SQA.QAdd) = left * right
liouvillian_product(left::Liouvillian, right::Liouvillian) = compose(left, right)
matrix_product(left::Matrix{ComplexF64}, right::Matrix{ComplexF64}) = left * right

function matrix_generator_norm(G::PeriodicGenerator{Matrix{ComplexF64}})
  isempty(keys(G)) && return 0.0
  return maximum(norm(G[harmonic]) for harmonic in keys(G))
end

function hamiltonian_wave_lhs(
  G::PeriodicGenerator{SQA.QAdd}, result::BlochReferenceResult, n::Int
)
  counts = BlochReferenceCounts()
  product_term = if iszero(n)
    G
  else
    reference_periodic_product(G, result.wave[n], hamiltonian_product, counts)
  end
  return SQA.simplify(product_term - im * derivative(result.wave[n + 1]))
end

function liouvillian_wave_lhs(
  G::PeriodicGenerator{Liouvillian}, result::BlochReferenceResult, n::Int
)
  counts = BlochReferenceCounts()
  product_term = if iszero(n)
    G
  else
    reference_periodic_product(G, result.wave[n], liouvillian_product, counts)
  end
  return SQA.simplify(product_term - derivative(result.wave[n + 1]))
end

function matrix_wave_lhs(
  G::PeriodicGenerator{Matrix{ComplexF64}}, result::BlochReferenceResult, n::Int
)
  counts = BlochReferenceCounts()
  product_term = if iszero(n)
    G
  else
    reference_periodic_product(G, result.wave[n], matrix_product, counts)
  end
  return product_term - im * derivative(result.wave[n + 1])
end

space = PauliSpace(:bloch_feshbach)
σx = Pauli(space, :σ, 1)
σy = Pauli(space, :σ, 2)
σz = Pauli(space, :σ, 3)
@variables w::Real

H = PeriodicGenerator(
  Dict(
    0 => 1 * σz, 1 => σx + im * σy, -1 => σx - im * σy, 2 => 2 * σx + σz, -2 => 2 * σx + σz
  ),
  w,
)

@testset "Bloch-Feshbach Hamiltonian reference recurrence" begin
  order = 3
  result = bloch_reference(
    H,
    order;
    product=hamiltonian_product,
    inverse_weight=harmonic -> 1 // harmonic,
    simplifier=SQA.simplify,
  )
  nonzero_harmonics = filter(!=(0), collect(keys(H)))

  @test length(result.effective) == order
  @test length(result.wave) == order - 1
  @test all(vanishes(time_average(X)) for X in result.wave)
  @test vanishes(result.effective[1] - H[0])

  for harmonic in nonzero_harmonics
    @test vanishes(result.wave[1][harmonic] - (1 // harmonic) * H[harmonic])
  end

  expected_B1 = zero(H[0])
  for harmonic in nonzero_harmonics
    expected_B1 += (1 // harmonic) * H[-harmonic] * H[harmonic]
  end
  @test vanishes(result.effective[2] - expected_B1)

  for harmonic in -4:4
    iszero(harmonic) && continue
    expected_X2 = zero(H[0])
    for inner_harmonic in nonzero_harmonics
      coefficient = 1 // (harmonic * inner_harmonic)
      expected_X2 += coefficient * H[harmonic - inner_harmonic] * H[inner_harmonic]
    end
    expected_X2 -= (1 // harmonic^2) * H[harmonic] * H[0]
    @test vanishes(result.wave[2][harmonic] - expected_X2)
  end

  expected_B2 = zero(H[0])
  for harmonic in nonzero_harmonics
    for inner_harmonic in nonzero_harmonics
      coefficient = 1 // (harmonic * inner_harmonic)
      expected_B2 +=
        coefficient * H[-harmonic] * H[harmonic - inner_harmonic] * H[inner_harmonic]
    end
    expected_B2 -= (1 // harmonic^2) * H[-harmonic] * H[harmonic] * H[0]
  end
  @test vanishes(result.effective[3] - expected_B2)

  for n in 0:(order - 2)
    lhs = hamiltonian_wave_lhs(H, result, n)
    rhs = reference_wave_equation_rhs(
      result, H, n, hamiltonian_product; simplifier=SQA.simplify
    )
    @test vanishes(lhs - rhs)
  end

  vv = floquet_expansion(H, VanVleck(), 2)
  @test vanishes(effective_component(vv, 1) - w^(-1) * result.effective[2])
  @test vanishes(micromotion(vv, 1) - (im * w^(-1)) * result.wave[1])
end

L = PeriodicGenerator(
  Dict(
    0 => hamiltonian_action(σz) + dissipator(σx),
    1 => hamiltonian_action(σx),
    -1 => dissipator(σy),
    2 => hamiltonian_action(σy),
    -2 => dissipator(σz),
  ),
  w,
)

@testset "Bloch-Feshbach Liouvillian reference recurrence" begin
  order = 3
  result = bloch_reference(
    L,
    order;
    product=liouvillian_product,
    inverse_weight=harmonic -> im // harmonic,
    simplifier=SQA.simplify,
  )
  nonzero_harmonics = filter(!=(0), collect(keys(L)))

  @test length(result.effective) == order
  @test length(result.wave) == order - 1
  @test all(liouvillian_vanishes(time_average(X)) for X in result.wave)
  @test liouvillian_vanishes(result.effective[1] - L[0])

  for harmonic in nonzero_harmonics
    expected = (im // harmonic) * L[harmonic]
    @test liouvillian_vanishes(result.wave[1][harmonic] - expected)
  end

  expected_B1 = zero(L[0])
  for harmonic in nonzero_harmonics
    expected_B1 += (im // harmonic) * compose(L[-harmonic], L[harmonic])
  end
  @test liouvillian_vanishes(result.effective[2] - expected_B1)

  for harmonic in -4:4
    iszero(harmonic) && continue
    expected_X2 = zero(L[0])
    for inner_harmonic in nonzero_harmonics
      coefficient = 1 // (harmonic * inner_harmonic)
      expected_X2 -=
        coefficient * compose(L[harmonic - inner_harmonic], L[inner_harmonic])
    end
    expected_X2 += (1 // harmonic^2) * compose(L[harmonic], L[0])
    @test liouvillian_vanishes(result.wave[2][harmonic] - expected_X2)
  end

  expected_B2 = zero(L[0])
  for harmonic in nonzero_harmonics
    for inner_harmonic in nonzero_harmonics
      coefficient = 1 // (harmonic * inner_harmonic)
      expected_B2 -= coefficient * compose(
        L[-harmonic], compose(L[harmonic - inner_harmonic], L[inner_harmonic])
      )
    end
    expected_B2 +=
      (1 // harmonic^2) * compose(L[-harmonic], compose(L[harmonic], L[0]))
  end
  @test liouvillian_vanishes(result.effective[3] - expected_B2)

  for n in 0:(order - 2)
    lhs = liouvillian_wave_lhs(L, result, n)
    rhs = reference_wave_equation_rhs(
      result, L, n, liouvillian_product; simplifier=SQA.simplify
    )
    @test liouvillian_vanishes(lhs - rhs)
  end

  vv = floquet_expansion(L, VanVleck(), 2)
  @test liouvillian_vanishes(effective_component(vv, 1) - w^(-1) * result.effective[2])
  @test liouvillian_vanishes(micromotion(vv, 1) - w^(-1) * result.wave[1])
end

M = PeriodicGenerator(
  Dict(
    0 => ComplexF64[0.4 1.0 + 0.2im; -0.3 + 0.1im -0.7],
    1 => ComplexF64[0.2 + 0.1im -0.6; 0.5im 0.8 - 0.2im],
    -1 => ComplexF64[-0.1 + 0.4im 0.3; 0.7 - 0.2im 0.5im],
    2 => ComplexF64[0.6 -0.2im; 0.9 + 0.1im -0.4],
    -2 => ComplexF64[-0.5im 0.8; -0.1 + 0.3im 0.2],
  ),
  w,
)

@testset "Bloch-Feshbach higher-order noncommuting matrix oracle" begin
  order = 6
  result = bloch_reference(
    M, order; product=matrix_product, inverse_weight=harmonic -> 1.0 / harmonic
  )
  lower = bloch_reference(
    M, 4; product=matrix_product, inverse_weight=harmonic -> 1.0 / harmonic
  )

  @test length(result.effective) == order
  @test length(result.wave) == order - 1
  @test all(norm(time_average(X)) <= 1.0e-14 for X in result.wave)

  for n in eachindex(lower.effective)
    @test norm(result.effective[n] - lower.effective[n]) <= 1.0e-12
  end
  for n in eachindex(lower.wave)
    for harmonic in union(keys(result.wave[n]), keys(lower.wave[n]))
      @test norm(result.wave[n][harmonic] - lower.wave[n][harmonic]) <= 1.0e-12
    end
  end

  for n in 0:(order - 2)
    lhs = matrix_wave_lhs(M, result, n)
    rhs = reference_wave_equation_rhs(result, M, n, matrix_product)
    @test matrix_generator_norm(lhs - rhs) <= 5.0e-11
  end
end

@testset "Bloch-Feshbach perturbative-order product count" begin
  for order in 1:8
    result = bloch_reference(
      zero(H),
      order;
      product=hamiltonian_product,
      inverse_weight=harmonic -> 1 // harmonic,
      simplifier=SQA.simplify,
    )
    @test result.counts.full_products == order - 1
    @test result.counts.static_products == order * (order - 1) ÷ 2
    @test result.counts.full_products + result.counts.static_products ==
      expected_bloch_series_products(order)
    @test iszero(result.counts.component_products)
  end
end

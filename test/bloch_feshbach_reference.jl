using Test
using FloquetExpansions
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
    H, order; product=hamiltonian_product, inverse_weight=harmonic -> 1 // harmonic
  )
  nonzero_harmonics = filter(!=(0), collect(keys(H)))

  @test length(result.effective) == order
  @test length(result.wave) == order - 1
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
    rhs = reference_wave_equation_rhs(result, H, n, hamiltonian_product)
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
    L, order; product=liouvillian_product, inverse_weight=harmonic -> im // harmonic
  )
  nonzero_harmonics = filter(!=(0), collect(keys(L)))

  @test length(result.effective) == order
  @test length(result.wave) == order - 1
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

  for n in 0:(order - 2)
    lhs = liouvillian_wave_lhs(L, result, n)
    rhs = reference_wave_equation_rhs(result, L, n, liouvillian_product)
    @test liouvillian_vanishes(lhs - rhs)
  end

  vv = floquet_expansion(L, VanVleck(), 2)
  @test liouvillian_vanishes(effective_component(vv, 1) - w^(-1) * result.effective[2])
  @test liouvillian_vanishes(micromotion(vv, 1) - w^(-1) * result.wave[1])
end

@testset "Bloch-Feshbach perturbative-order product count" begin
  for order in 1:8
    result = bloch_reference(
      zero(H), order; product=hamiltonian_product, inverse_weight=harmonic -> 1 // harmonic
    )
    @test result.counts.full_products == order - 1
    @test result.counts.static_products == order * (order - 1) ÷ 2
    @test result.counts.full_products + result.counts.static_products ==
      expected_bloch_series_products(order)
    @test iszero(result.counts.component_products)
  end
end

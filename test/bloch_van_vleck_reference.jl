using Test
using FloquetExpansions
using LinearAlgebra: norm
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

include(joinpath(@__DIR__, "helpers", "shared.jl"))
include(joinpath(@__DIR__, "helpers", "bloch_van_vleck_reference.jl"))

hamiltonian_product_vv(left::SQA.QAdd, right::SQA.QAdd) = left * right
liouvillian_product_vv(left::Liouvillian, right::Liouvillian) = compose(left, right)
matrix_product_vv(left::Matrix{ComplexF64}, right::Matrix{ComplexF64}) = left * right

liouvillian_vanishes_vv(L::Liouvillian) = iszero(SQA.simplify(L))
function liouvillian_vanishes_vv(G::PeriodicGenerator{Liouvillian})
  return all(liouvillian_vanishes_vv(G[harmonic]) for harmonic in keys(G))
end

function matrix_generator_norm_vv(G::PeriodicGenerator{Matrix{ComplexF64}})
  isempty(keys(G)) && return 0.0
  return maximum(norm(G[harmonic]) for harmonic in keys(G))
end

space = NLevelSpace(:bloch_van_vleck, 3)
σ11 = Transition(space, :σ, 1, 1)
σ22 = Transition(space, :σ, 2, 2)
σ33 = Transition(space, :σ, 3, 3)
σ12 = Transition(space, :σ, 1, 2)
σ23 = Transition(space, :σ, 2, 3)
σ31 = Transition(space, :σ, 3, 1)
σ13 = Transition(space, :σ, 1, 3)
σ21 = Transition(space, :σ, 2, 1)
σ32 = Transition(space, :σ, 3, 2)
@variables w_vv::Real

H0_vv = 2 * σ11 - σ22 + (3 // 2) * σ33 + σ12 + σ12'
H1_vv = σ12 + 2 * σ23 + im * σ31
H2_vv = 2 * σ13 - σ21 + im * σ32
H_vv = PeriodicGenerator(
  Dict(0 => H0_vv, 1 => H1_vv, -1 => H1_vv', 2 => H2_vv, -2 => H2_vv'), w_vv
)

@testset "Bloch-Feshbach converts to canonical Hamiltonian Van Vleck" begin
  order = 3
  bloch = bloch_reference(
    H_vv,
    order;
    product=hamiltonian_product_vv,
    inverse_weight=harmonic -> 1 // harmonic,
    simplifier=SQA.simplify,
  )
  converted = bloch_van_vleck_reference(
    bloch, H_vv; product=hamiltonian_product_vv, simplifier=SQA.simplify
  )
  vv = floquet_expansion(H_vv, VanVleck(), order)

  @test length(converted.log_embedding) == order - 1
  @test all(vanishes(time_average(Gn)) for Gn in converted.log_embedding)
  @test converted.counts.log_products == expected_mercator_products(order - 1)

  for n in 1:(order - 1)
    @test vanishes(converted.effective[n + 1] - w_vv^n * effective_component(vv, n))
    kick_n = im * converted.log_embedding[n]
    @test vanishes(kick_n - w_vv^n * micromotion(vv, n))
  end
end

L_vv = PeriodicGenerator(
  Dict(
    0 => hamiltonian_action(σ11 - σ33) + dissipator(σ12 + σ23),
    1 => hamiltonian_action(σ12 + σ21),
    -1 => dissipator(σ23 + σ31),
    2 => hamiltonian_action(σ13 + σ31),
    -2 => dissipator(σ12 - im * σ23),
  ),
  w_vv,
)

@testset "Bloch-Feshbach converts to canonical Liouvillian Van Vleck" begin
  order = 3
  bloch = bloch_reference(
    L_vv,
    order;
    product=liouvillian_product_vv,
    inverse_weight=harmonic -> im // harmonic,
    simplifier=SQA.simplify,
  )
  converted = bloch_van_vleck_reference(
    bloch, L_vv; product=liouvillian_product_vv, simplifier=SQA.simplify
  )
  vv = floquet_expansion(L_vv, VanVleck(), order)

  @test length(converted.log_embedding) == order - 1
  @test all(liouvillian_vanishes_vv(time_average(Gn)) for Gn in converted.log_embedding)
  @test converted.counts.log_products == expected_mercator_products(order - 1)

  for n in 1:(order - 1)
    @test liouvillian_vanishes_vv(
      converted.effective[n + 1] - w_vv^n * effective_component(vv, n)
    )
    @test liouvillian_vanishes_vv(converted.log_embedding[n] - w_vv^n * micromotion(vv, n))
  end
end

M_vv = PeriodicGenerator(
  Dict(
    0 => ComplexF64[0.3 0.7 - 0.1im; -0.2 + 0.4im -0.5],
    1 => ComplexF64[0.4im 0.8; -0.3 0.2 + 0.1im],
    -1 => ComplexF64[0.6 -0.2im; 0.5im -0.1],
  ),
  w_vv,
)

@testset "Dense noncommuting Van Vleck conversion through order six" begin
  order = 6
  tolerance = 5.0e-9

  bloch_hamiltonian = bloch_reference(
    M_vv, order; product=matrix_product_vv, inverse_weight=harmonic -> 1.0 / harmonic
  )
  converted_hamiltonian = bloch_van_vleck_reference(
    bloch_hamiltonian, M_vv; product=matrix_product_vv
  )
  hori_hamiltonian = hori_deprit_reference(
    M_vv, order; product=matrix_product_vv, phase=im
  )

  @test converted_hamiltonian.counts.log_products == expected_mercator_products(order - 1)
  @test converted_hamiltonian.counts.log_products == 20
  @test all(norm(time_average(Gn)) <= 1.0e-12 for Gn in converted_hamiltonian.log_embedding)

  for n in 1:order
    @test norm(converted_hamiltonian.effective[n] - hori_hamiltonian.effective[n]) <= tolerance
  end
  for n in 1:(order - 1)
    kick_residual = im * converted_hamiltonian.log_embedding[n] - hori_hamiltonian.kick[n]
    @test matrix_generator_norm_vv(kick_residual) <= tolerance
  end

  bloch_map = bloch_reference(
    M_vv, order; product=matrix_product_vv, inverse_weight=harmonic -> im / harmonic
  )
  converted_map = bloch_van_vleck_reference(bloch_map, M_vv; product=matrix_product_vv)
  hori_map = hori_deprit_reference(M_vv, order; product=matrix_product_vv, phase=-1)

  @test all(norm(time_average(Gn)) <= 1.0e-12 for Gn in converted_map.log_embedding)
  for n in 1:order
    @test norm(converted_map.effective[n] - hori_map.effective[n]) <= tolerance
  end
  for n in 1:(order - 1)
    @test matrix_generator_norm_vv(converted_map.log_embedding[n] - hori_map.kick[n]) <= tolerance
  end

  lower = bloch_van_vleck_reference(
    bloch_reference(
      M_vv, 4; product=matrix_product_vv, inverse_weight=harmonic -> 1.0 / harmonic
    ),
    M_vv;
    product=matrix_product_vv,
  )
  for n in eachindex(lower.effective)
    @test norm(lower.effective[n] - converted_hamiltonian.effective[n]) <= tolerance
  end
  for n in eachindex(lower.log_embedding)
    @test matrix_generator_norm_vv(
      lower.log_embedding[n] - converted_hamiltonian.log_embedding[n]
    ) <= tolerance
  end
end

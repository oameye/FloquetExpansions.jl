using Test
using Symbolics: @variables
using FloquetExpansions
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions

@variables ω::Real t::Real Δ::Real χ::Real κ::Real

function driven_kerr(dissipative::Bool)
  h = FockSpace(:c)
  a = Destroy(h, :a)
  H = Δ * a' * a + χ * a' * a * a' * a + (2 // 5) * (a + a') * cos(ω * t)
  L = dissipative ? liouvillian(H; channels=(jump(a, κ),)) : liouvillian(H)
  return harmonics(L, ω, t), a
end

const KERR_VALUES = Dict(Δ => 1 // 2, χ => 3 // 10, κ => 4 // 5)

# The Van Vleck reference is computed in Float64, so its lowering is compared with a tolerance.
function nearly_same(A::FE.AlgebraSuperoperator, B::FE.AlgebraSuperoperator)
  return all(abs(ComplexF64(c)) <= 1e-12 for c in values((A + (-1) * B).terms))
end

@testset "native expansion equals Van Vleck wherever no static gauge is needed" begin
  for dissipative in (false, true), algorithm in (BlochFeshbach(), HoriDeprit())
    G, _ = driven_kerr(dissipative)
    data = FE.native_expansion_data(algorithm, G, 3, KERR_VALUES, 3)
    reference = floquet_expansion(G, VanVleck(; algorithm), 4)
    @test all(iszero, data.recurrence.S)
    for n in 0:3
      @test nearly_same(
        data.recurrence.E[n + 1],
        FE.lower_liouvillian(data.lowering, reference.effective_components[n + 1]),
      )
    end
  end
end

@testset "cutoff-free driven Kerr births one exact channel at order four" begin
  G, a = driven_kerr(true)
  bf = FE.native_expansion_data(BlochFeshbach(), G, 4, KERR_VALUES, 3)
  hd = FE.native_expansion_data(HoriDeprit(), G, 4, KERR_VALUES, 3)
  @test bf.recurrence.E == hd.recurrence.E
  newborn = bf.recurrence.channels[2:end]
  @test length(newborn) == 1
  @test newborn[1].onset == 4
  @test newborn[1].weight == 72 // 3125
  jump_operator = FE.frame_operator(bf.representation, newborn[1].coefficients[1])
  @test jump_operator == FE.lower_qadd(bf.lowering, FE.qadd(a * a - 2 * a' * a))
end

@testset "SQA lowering is exact across Fock, level and Pauli sites" begin
  h = FockSpace(:c) ⊗ NLevelSpace(:q, 2) ⊗ PauliSpace(:p)
  a = Destroy(h, :a)
  σ = Transition(h, :σ, 1, 2)
  sx = Pauli(h, :s, 1)
  sz = Pauli(h, :s, 3)
  H = Δ * a' * a + (1 // 2) * sz + (a' * σ + a * σ') + (2 // 5) * (a + a') * cos(ω * t) * sx
  G = harmonics(liouvillian(H; channels=(collapse(sqrt(κ) * a),)), ω, t)
  lowering = FE.SQALowering{Complex{Rational{Int128}}}(G, KERR_VALUES)
  identity = FE.algebra_identity(lowering.algebra)
  @test FE.lower_qadd(lowering, FE.qadd(a * a' - a' * a)) == identity
  @test FE.lower_qadd(lowering, FE.qadd(sx * sx)) == identity
  X = FE.lower_qadd(lowering, FE.qadd(a' * σ + sz))
  @test FE.lower_qadd(lowering, FE.lift_operator(lowering, X)) == X
  incomplete = FE.SQALowering{Complex{Rational{Int128}}}(G, Dict(Δ => 1 // 2))
  @test_throws ArgumentError FE.lower_generator(incomplete, G)
end

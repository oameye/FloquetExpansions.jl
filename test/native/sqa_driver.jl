using Test
using Symbolics: @variables
using FloquetExpansions
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions

@variables ω::Real t::Real

function driven_kerr(dissipative::Bool)
  h = FockSpace(:c)
  a = Destroy(h, :a)
  H = (1 // 2) * a' * a + (3 // 10) * a' * a * a' * a + (2 // 5) * (a + a') * cos(ω * t)
  L = hamiltonian_action(H)
  dissipative && (L += (4 // 5) * dissipator(a))
  return harmonics(L, ω, t), a
end

@testset "native expansion equals Van Vleck wherever no static gauge is needed" begin
  for dissipative in (false, true), algorithm in (BlochFeshbach(), HoriDeprit())
    G, _ = driven_kerr(dissipative)
    data = FE.native_expansion_data(algorithm, G, 3, 3)
    reference = floquet_expansion(G, VanVleck(; algorithm), 4)
    @test all(iszero, data.recurrence.S)
    for n in 0:3
      @test data.recurrence.E[n + 1] ==
        FE.lower_liouvillian(data.lowering, reference.effective_components[n + 1])
    end
  end
end

@testset "cutoff-free driven Kerr births one exact channel at order four" begin
  G, a = driven_kerr(true)
  bf = FE.native_expansion_data(BlochFeshbach(), G, 4, 3)
  hd = FE.native_expansion_data(HoriDeprit(), G, 4, 3)
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
  H =
    (1 // 2) * a' * a +
    (1 // 2) * sz +
    (a' * σ + a * σ') +
    (2 // 5) * (a + a') * cos(ω * t) * sx
  G = harmonics(hamiltonian_action(H) + dissipator((2 // 3) * a), ω, t)
  lowering = FE.SQALowering{Complex{Rational{Int128}}}(G)
  identity = FE.algebra_identity(lowering.algebra)
  @test FE.lower_qadd(lowering, FE.qadd(a * a' - a' * a)) == identity
  @test FE.lower_qadd(lowering, FE.qadd(sx * sx)) == identity
  X = FE.lower_qadd(lowering, FE.qadd(a' * σ + sz))
  @test FE.lower_qadd(lowering, FE.lift_operator(lowering, X)) == X
end

@testset "the exact driver falls back to BigInt scalars when Int128 overflows" begin
  h = FockSpace(:c)
  a = Destroy(h, :a)
  H =
    (1 // 1000000) * a' * a +
    (3 // 1000007) * a' * a * a' * a +
    (2 // 1000009) * (a + a') * cos(ω * t)
  G = harmonics(hamiltonian_action(H) + (4 // 1000011) * dissipator(a), ω, t)
  @test_throws OverflowError FE.native_expansion_data(
    Complex{Rational{Int128}}, HoriDeprit(), G, 2, 3
  )
  data = FE.native_expansion_data(HoriDeprit(), G, 2, 3)
  exact = FE.native_expansion_data(Complex{Rational{BigInt}}, HoriDeprit(), G, 2, 3)
  @test data isa FE.NativeExpansionData{Complex{Rational{BigInt}},Rational{BigInt}}
  @test data.recurrence.E == exact.recurrence.E
end

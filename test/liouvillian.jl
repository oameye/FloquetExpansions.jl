using Test
using FloquetExpansions
import SecondQuantizedAlgebra as SQA
using Symbolics: @variables

space = FockSpace(:liouvillian)
a = Destroy(space, :a)
@variables γ::Real

@testset "Liouvillian channel constructors" begin
  H = a' * a
  coherent = hamiltonian_action(H)
  complete = dissipator(a)
  weighted = γ * dissipator(a')

  L = Liouvillian(H; collapse_operators=(a,), jumps=(a',), rates=(γ,))

  @test L == coherent + complete + weighted
  @test L isa Liouvillian
  @test @inferred(hamiltonian_action(H)) isa Liouvillian
  @test @inferred(dissipator(a)) isa Liouvillian
  @test @inferred(γ * dissipator(a')) isa Liouvillian
  @test @inferred(Liouvillian(H; jumps=(a,), rates=(γ,))) isa Liouvillian
end

@testset "Liouvillian arithmetic collects equal actions" begin
  L = hamiltonian_action(a' * a)

  @test iszero(L - L)
  @test zero(L) + L == L
  @test 2 * L == L + L
  @test L * 2 == L + L
  @test SQA.simplify(L - L) == zero(L)
end

@testset "Liouvillian composition is map composition" begin
  coherent = hamiltonian_action(a' * a)
  dissipative = dissipator(a)
  composed = compose(coherent, dissipative)

  @test composed isa Liouvillian
  @test !iszero(composed)
  @test compose(coherent, zero(dissipative)) == zero(coherent)
  @test compose(zero(coherent), dissipative) == zero(coherent)
  @test SQA.commutator(coherent, dissipative) ==
    compose(coherent, dissipative) - compose(dissipative, coherent)
  @test compose(coherent, dissipative) != compose(dissipative, coherent)
end

@testset "Liouvillian constructor validates paired channels" begin
  @test_throws ArgumentError Liouvillian(a; jumps=(a,), rates=())
  @test_throws ArgumentError Liouvillian(a; jumps=(a,), rates=(1, 2))
end

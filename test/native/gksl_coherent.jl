using Test
using Symbolics: @variables
using FloquetExpansions

@variables ω::Real t::Real

@testset "GKSLNormalForm accepts coherent Liouvillian input" begin
  a = Destroy(FockSpace(:c), :a)
  for H in (a' * a, a' * a + (a + a') * cos(ω * t))
    G = harmonics(hamiltonian_action(H), ω, t)
    for algorithm in (HoriDeprit(), BlochFeshbach())
      native = floquet_expansion(G, GKSLNormalForm(; algorithm), 3)
      reference = floquet_expansion(G, VanVleck(; algorithm), 3)
      @test isempty(channels(native))
      @test size(kossakowski(native, dissipative_frame(native))) == (0, 0)
      @test size(kossakowski(native)) == (0, 0)
      @test iszero(simplify(hamiltonian(native) - hamiltonian(reference)))
      @test iszero(simplify(effective_generator(native) - effective_generator(reference)))
      @test effective_generator(native) == hamiltonian_action(hamiltonian(native))
      for n in 0:2
        @test size(kossakowski_component(native, n)) == (0, 0)
        @test iszero(
          simplify(effective_component(native, n) - effective_component(reference, n))
        )
      end
      for n in 1:2
        @test iszero(simplify(micromotion(native, n) - micromotion(reference, n)))
      end
    end
  end
end

@testset "GKSLNormalForm accepts zero Liouvillian input" begin
  for G in (
    PeriodicGenerator(Dict(0 => zero(Liouvillian)), ω),
    PeriodicGenerator(Dict{Int,Liouvillian}(), ω),
  )
    for algorithm in (HoriDeprit(), BlochFeshbach()), N in (1, 2, 3)
      native = floquet_expansion(G, GKSLNormalForm(; algorithm), N)
      @test order(native) == N
      @test iszero(effective_generator(native))
      @test iszero(hamiltonian(native))
      @test iszero(micromotion(native))
      @test isempty(channels(native))
      @test size(kossakowski(native, dissipative_frame(native))) == (0, 0)
      @test size(kossakowski(native)) == (0, 0)
      for n in 0:(N - 1)
        @test iszero(effective_component(native, n))
        @test size(kossakowski_component(native, n)) == (0, 0)
      end
    end
  end
end

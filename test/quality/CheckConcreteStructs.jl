using Test
using FloquetExpansions
using CheckConcreteStructs: all_concrete

@testset "no abstract struct fields" begin
  @test all_concrete(FloquetExpansions)
end

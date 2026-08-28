using Test
using FloquetExpansions
using CheckConcreteStructs: all_concrete

# CLAUDE.md: no abstract-typed fields anywhere. `PeriodicOperator.components` is
# `Dict{Int,QAdd}` and `FloquetExpansion`'s vectors are concrete for the same reason.
@testset "no abstract struct fields" begin
  @test all_concrete(FloquetExpansions)
end

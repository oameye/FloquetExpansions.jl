using Test
using FloquetExpansions
using JET: JET

@testset "JET" begin
  JET.test_package(FloquetExpansions; target_defined_modules=true)
end

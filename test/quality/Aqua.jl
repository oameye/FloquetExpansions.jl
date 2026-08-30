using Test
using FloquetExpansions
using Aqua: Aqua

@testset "Aqua" begin
  Aqua.test_all(FloquetExpansions; project_extras=false)
end

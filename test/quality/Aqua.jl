using Test
using FloquetExpansions
using Aqua: Aqua

@testset "Aqua" begin
  # `project_extras` is for the legacy `[extras]` layout; this package uses test/Project.toml.
  Aqua.test_all(FloquetExpansions; project_extras=false)
end

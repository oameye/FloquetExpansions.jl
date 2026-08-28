using Test
using FloquetExpansions
using ExplicitImports:
  ExplicitImports,
  check_all_explicit_imports_via_owners,
  check_all_qualified_accesses_via_owners,
  check_no_implicit_imports,
  check_no_self_qualified_accesses,
  check_no_stale_explicit_imports

@testset "ExplicitImports" begin
  # `@reexport using SecondQuantizedAlgebra` is deliberate: the operator algebra is the user's
  # working vocabulary and re-exporting it is the point. It is the one implicit import allowed.
  @test check_no_implicit_imports(
    FloquetExpansions; skip=(Base, Core, SecondQuantizedAlgebra)
  ) === nothing
  @test check_no_stale_explicit_imports(FloquetExpansions) === nothing
  @test check_all_explicit_imports_via_owners(FloquetExpansions) === nothing
  @test check_all_qualified_accesses_via_owners(FloquetExpansions) === nothing
  @test check_no_self_qualified_accesses(FloquetExpansions) === nothing
end

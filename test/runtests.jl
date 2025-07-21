using Test
using FloquetExpansions

if isempty(VERSION.prerelease)
    @testset "Code linting" begin
        using JET
        JET.test_package(FloquetExpansions; target_defined_modules=true)
        # rep = report_package("FloquetExpansions")
        # @show rep
        # @test length(JET.get_reports(rep)) <= 2
        # @test_broken length(JET.get_reports(rep)) == 0
    end
end

@testset "ExplicitImports" begin
    using ExplicitImports
    @test check_no_implicit_imports(FloquetExpansions) == nothing
    @test check_all_explicit_imports_via_owners(FloquetExpansions) == nothing
    # @test check_all_explicit_imports_are_public(FloquetExpansions) == nothing
    @test check_no_stale_explicit_imports(FloquetExpansions) == nothing
    @test check_all_qualified_accesses_via_owners(FloquetExpansions) == nothing
    # @test check_all_qualified_accesses_are_public(FloquetExpansions) == nothing
    @test check_no_self_qualified_accesses(FloquetExpansions) == nothing
end

@testset "best practices" begin
    using Aqua

    Aqua.test_ambiguities([FloquetExpansions]; broken=false)
    Aqua.test_all(FloquetExpansions; ambiguities=false)
end

@testset "Concretely typed" begin
    import FloquetExpansions as FE
    using CheckConcreteStructs
end

# Include all test files

@testset "rotation" begin
    include("utils.jl")
end

@testset "rotation" begin
    include("rotate.jl")
end
@testset "Kerr Resonator" begin
    include("kerr_resonator.jl")
end

@testset "Documentation" begin
    using Documenter
    DocMeta.setdocmeta!(
        FloquetExpansions, :DocTestSetup, :(using FloquetExpansions); recursive=true
    )
    Documenter.doctest(FloquetExpansions)
end

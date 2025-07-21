using FloquetExpansions
using SymbolicUtils
using Test

@testset "remove constants" begin
    # Hilbert space
    h = FockSpace(:cavity)

    @syms Δ F κ ω ω₀ α
    @qnumbers a::Destroy(h)
    @syms t::Real

    Ht = ω₀ * a' * a + α * (a' + a)^4 / 4 + F * (a' + a) * cos(ω * t)
    rot = rotate(Ht, a, ω, t)
    rot_ = FloquetExpansions.remove_constants!(deepcopy(rot))
    @test length(rot_.arguments) == 14
    @test length(rot_.arguments) !== length(rot.arguments)
end

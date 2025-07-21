using FloquetExpansions
using SymbolicUtils
using Test
using QuestBase: @eqtest

@testset "external driving" begin
    # Hilbert space
    h = FockSpace(:cavity)

    @syms Δ F κ ω ω₀ α
    @qnumbers a::Destroy(h)
    @syms t::Real

    Ht = (F * (a' + a) * cos(ω * t))
    rot = rotate(Ht, a, ω, t)

    QE_operator = quasienergy_operator(rot, ω, t)

    @test Set(keys(QE_operator)) == Set([0, 2, -2])
    @eqtest QE_operator[2].arguments[1] == (1//2 * F * a)
end

@testset "Kerr resonator" begin
    # Hilbert space
    h = FockSpace(:cavity)

    @syms Δ F κ ω ω₀ α
    @qnumbers a::Destroy(h)
    @syms t::Real

    #
    Ht = simplify(ω₀ * a' * a + α * (a' + a)^4 / 4 + F * (a' + a) * cos(ω * t))
    rot = rotate(Ht, a, ω, t)
    QE_operator = quasienergy_operator(rot, ω, t)

    @test Set(keys(QE_operator)) == Set([0, 2, -2, 4, -4])

    @eqtest simplify(QE_operator[0] - (
        ω₀ * (a' * a) +
        3.0α * (a' * a) +
        0.5F * (a') +
        0.5F * (a) +
        1.5α * (a' * a' * a * a) +
        -ω * (a' * a)
    )) == 0
    @eqtest QE_operator[4].arguments[1] == (0.25α * (a * a * a * a))
end

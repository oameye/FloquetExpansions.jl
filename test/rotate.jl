using FloquetExpansions
using SymbolicUtils
using QuestBase: @eqtest
using Test

# Hilbert space
h = FockSpace(:cavity)

@cnumbers Δ F κ ω ω₀
@qnumbers a::Destroy(h)
@syms t::Real

@testset begin
    expr = simplify((a' + a))
    rot = rotate(expr, a, ω, t)

    # @test simplify(rot.arguments[end] + ω * a' * a) == 0
    @test iszero(simplify(rot - (exp((0 - 1im)*t*ω)*(a')+exp((0 + 1im)*t*ω)*(a)+-ω*(a'*a))))
    # @eqtest rot ==
    # (exp((0 - 1im) * t * ω) * (a') + exp((0 + 1im) * t * ω) * (a) - ω * (a' * a))
end

@testset begin
    expr = simplify(F * (a' + a))
    rot = rotate(expr, a, ω, t)

    # @test simplify(rot.arguments[end] + ω * a' * a) == 0
    # @test string(rot) == "(F*exp((0 - 1im)*t*ω)*(a′)+F*exp((0 + 1im)*t*ω)*(a)+-ω*(a′*a))"
    @test iszero(
        simplify(
            rot - (
                F * exp((0 - 1im) * t * ω) * (a') + F * exp((0 + 1im) * t * ω) * (a) -
                ω * (a' * a)
            ),
        ),
    )
end

@testset begin
    using SecondQuantizedAlgebra: average, undo_average
    expr = simplify(ω₀ * (a' * a))
    rot = simplify(rotate(expr, a, ω, t))

    @test simplify(average(rot) - average(-ω * (a' * a) + ω₀ * (a' * a))) == 0
    # @test string(rot) == "(-ω*(a′*a)+ω₀*(a′*a))"
    # @eqtest rot == (-ω*(a'*a) + ω₀*(a'*a))
end

@testset begin
    expr = simplify((a' + a)^4)
    rot = rotate(expr, a, ω, t)

    term = average(
        6 * (a' * a' * a * a) +
        3 +
        12 * (a' * a) +
        exp((0 + 4im) * t * ω) * (a * a * a * a) +
        6exp((0 + 2im) * t * ω) * (a * a) +
        4exp((0 + 2im) * t * ω) * (a' * a * a * a) +
        6exp((0 - 2im) * t * ω) * (a' * a') +
        4exp((0 - 2im) * t * ω) * (a' * a' * a' * a) +
        exp((0 - 4im) * t * ω) * (a' * a' * a' * a') - ω * (a' * a),
    )
    @test simplify(average(rot) - term) == 0
end

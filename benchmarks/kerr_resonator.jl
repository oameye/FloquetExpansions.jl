function benchmark_kerr_resonator!(SUITE)
    using FloquetExpansions

    h = FockSpace(:cavity)
    @rnumbers Δ F κ ω ω₀ α
    @qnumbers a::Destroy(h)

    Ht = simplify(ω₀ * a' * a + α * (a' + a)^4 / 4 + F * (a' + a) * cos(ω * t))
    rotating_wave_approximation(Ht, a, ω, t)

    return SUITE["Operation"]["Rotate"]["Kerr resonator"] = @benchmarkable rotate(
        $Ht, $a, $ω, $t
    ) seconds = 10

    rot = rotate(Ht, a, ω, t)
    return SUITE["Construction"]["Quasienergy operator"]["Kerr resonator"] = @benchmarkable quasienergy_operator(
        $rot, $ω, $t
    ) seconds = 10

    return SUITE["Floquet Expansion"]["Rotating Wave Approximation"]["Kerr resonator"] = @benchmarkable rotating_wave_approximation(
        $Ht, $a, $ω, $t
    ) seconds = 10
end

function rotating_wave_approximation(
    H::QTerm, a::QSym, ω=SQA.cnumber(:ω), t=SQA.rnumber(:t)
)
    rot = rotate(H, a, ω, t)
    components = get_fourier_components(rot, ω, t)
    return components[0]
end

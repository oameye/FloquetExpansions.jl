function rotating_wave_approximation(
    H::QTerm, a::QSym, ω=SQA.cnumber(:ω), t=SQA.rnumber(:t)
)
    Hrot = rotate(H, a, ω, t)
    Hrot_av = remove_constants(SQA.average(Hrot))
    components = get_fourier_components(Hrot_av, ω, t)
    return components[0]
end

function rotating_wave_approximation(H::QTerm, a::QSym, ω, t)
    components = quasienergy_operator(H, a, ω, t)
    return components[0]
end

function quasienergy_operator(H::QTerm, a::QSym, ω, t)
    rot = rotate(H, a, ω, t)
    components = get_fourier_components(rot, ω, t)
    return components
end

function quasienergy_operator(H::QTerm, ω, t)
    components = get_fourier_components(H, ω, t)
    return components
end

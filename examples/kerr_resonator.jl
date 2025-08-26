# using FloquetExpansions

# h = FockSpace(:cavity)
# @rnumbers Δ F κ ω ω₀ α t
# @qnumbers a::Destroy(h)

# Ht = ω₀ * a' * a + α * (a' + a)^4 / 4 + F * (a' + a) * cos(ω * t)

# rotating_wave_approximation(Ht, a, ω, t)

# #

# Q = quasienergy_operator(Ht, a, ω, t)


using VanVleckRecursion

H = Terms([
    Term(; rotating=0),    # H₀ (static part)
    Term(; rotating=1),     # H₁cos(ωt) (oscillating part)
    # Term(; rotating=2),     # H₁cos(ωt) (oscillating part)
])

set_hamiltonian!(H)

# K(2)
# term= K(2)
# latex(term)

term= K(3).terms[end]
latex(term)

using FloquetExpansions
using SymbolicUtils
using Test

# Hilbert space
h = FockSpace(:cavity)

@syms Δ F κ ω ω₀ α
@qnumbers a::Destroy(h)
@syms t::Real

#
Ht = simplify(ω₀ * a' * a + α * (a' + a)^4 / 4 + F * (a' + a) * cos(ω * t))
rot = rotate(Ht, a, ω, t)
rot = FloquetExpansions.remove_constants!(deepvopy(rot))

get_fourier_components(rot, ω, t)

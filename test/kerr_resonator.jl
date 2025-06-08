using FloquetExpansions
using SymbolicUtils

# Hilbert space
h = FockSpace(:cavity)

@syms Δ F κ ω ω₀ α
@qnumbers a::Destroy(h)
@syms t::Real

#

Ht = simplify(ω₀ * a' * a + α * (a' + a)^4 / 4 + F * (a' + a) * cos(ω * t))
rot = rotate(Ht, a, ω, t)
rot_av = FloquetExpansions.remove_constants(FloquetExpansions.average(rot))
term = arguments(rot_av)[1]
@test SymbolicUtils.ismul(term)

```@meta
EditURL = "../../../examples/kerr_resonator.jl"
```

````@example kerr_resonator
using FloquetExpansions

h = FockSpace(:cavity)
@rnumbers Δ F κ ω ω₀ α t
@qnumbers a::Destroy(h)

Ht = ω₀ * a' * a + α * (a' + a)^4 / 4 + F * (a' + a) * cos(ω * t)

rotating_wave_approximation(Ht, a, ω, t)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


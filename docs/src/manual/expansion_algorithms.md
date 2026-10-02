```@meta
CurrentModule = FloquetExpansions
```

# [Expansion algorithms](@id expansion-algorithms-manual)

After the gauge, the second choice in a Floquet expansion is the algorithm that computes it.

```@docs
ExpansionAlgorithm
```

## Choosing an algorithm

```@docs
HoriDeprit
BlochFeshbach
```

This example benchmarks prepared qubit generators. See the
[Bloch/Feshbach theory page](@ref bloch-feshbach-theory) for the recurrence-scaling analysis.

```@example algorithm-cost
using BenchmarkTools: @benchmark, median
using FloquetExpansions
using Plots
using Symbolics: @variables

space = PauliSpace(:algorithm_cost)
σx = Pauli(space, :σ, 1)
σz = Pauli(space, :σ, 3)
@variables ω::Real t::Real Δ::Real A::Real γ::Real
H = (Δ / 2) * σz + A * cos(ω * t) * σx
periodic_hamiltonian = harmonics(H, ω, t)
periodic_liouvillian = harmonics(liouvillian(H; channels=(jump(σx, γ),)), ω, t)

function median_runtime_ms(generator, algorithm, order)
  gauge = VanVleck(; algorithm)
  trial = @benchmark floquet_expansion($generator, $gauge, $order) samples=7 evals=1
  return median(trial).time / 1e6
end

hamiltonian_orders = 1:10
liouvillian_orders = 1:6
hamiltonian_times = hcat(
  [median_runtime_ms(periodic_hamiltonian, HoriDeprit(), n) for n in hamiltonian_orders],
  [median_runtime_ms(periodic_hamiltonian, BlochFeshbach(), n) for n in hamiltonian_orders],
)
liouvillian_times = hcat(
  [median_runtime_ms(periodic_liouvillian, HoriDeprit(), n) for n in liouvillian_orders],
  [median_runtime_ms(periodic_liouvillian, BlochFeshbach(), n) for n in liouvillian_orders],
)

hamiltonian_plot = plot(
  hamiltonian_orders, hamiltonian_times;
  xscale=:log10, yscale=:log10, xticks=(hamiltonian_orders, string.(hamiltonian_orders)),
  marker=:circle, label=["Hori–Deprit" "Bloch/Feshbach"], title="Hamiltonian",
  xlabel="order N", ylabel="median runtime (ms)", legend=:topleft,
)
liouvillian_plot = plot(
  liouvillian_orders, liouvillian_times;
  xscale=:log10, yscale=:log10, xticks=(liouvillian_orders, string.(liouvillian_orders)),
  marker=:circle, label=["Hori–Deprit" "Bloch/Feshbach"], title="Liouvillian",
  xlabel="order N", ylabel="median runtime (ms)", legend=:topleft,
)
plot(hamiltonian_plot, liouvillian_plot; layout=(1, 2), size=(1100, 500))
```

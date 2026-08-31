# # Counting pump photons: a mean-field comparison

using FloquetExpansions

h_a = FockSpace(:cavity)

@qnumbers a::Destroy(h_a)

@variables ω₀::Real ω::Real α::Real F::Real γ::Real
@variables t::Real

x_a = (a + a') / sqrt(2 * ω₀)
p_a = im * sqrt(ω₀ / 2) * (a' - a)
function duffing_hamiltonian(x, p)
  return p^2 // 2 + ω₀^2 * x^2 // 2 + α * x^4 // 4 - F * x * exponential_form(cos(ω * t))
end

duffing = duffing_hamiltonian(x_a, p_a)

#

rot_a = Rotation(a, ω * t, t)

H_a = simplify(transform(duffing, rot_a))

#

Q_a = harmonics(H_a, ω, t)

#

rwa_a = floquet_expansion(Q_a, VanVleck(), 1)
H_eff_a = simplify(effective_hamiltonian(rwa_a))

#

h_b = FockSpace(:cavity)
@qnumbers b::Destroy(h_b)

x_b = (b + b') / sqrt(2 * ω)
p_b = im * sqrt(ω / 2) * (b' - b)

rot_b = Rotation(b, ω * t, t)

H_b = simplify(transform(duffing_hamiltonian(x_b, p_b), rot_b))

#

Q_b = harmonics(H_b, ω, t)
rwa_b = floquet_expansion(Q_b, VanVleck(), 1)
H_eff_b = simplify(effective_hamiltonian(rwa_b))

#

using QuantumCumulants, ModelingToolkitBase

eqs_a = meanfield([a], H_eff_a, [a]; rates=[γ], order=1, iv=t)
eqs_b = meanfield([b], H_eff_b, [b]; rates=[γ], order=1, iv=t)

#

sys_a = mtkcompile(System(eqs_a; name=:a_order_1))
sys_b = mtkcompile(System(eqs_b; name=:b_order_1))

u0_a = initial_values(eqs_a; defaults=Dict(average(a) => 0.0 + 0.0im))
u0_b = initial_values(eqs_b; defaults=Dict(average(b) => 0.0 + 0.0im))

#
using OrdinaryDiffEq: Tsit5, solve

params = Dict(ω₀ => 1.0, ω => 2, α => 1.0, F => 0.01, γ => 0.005)
tspan = (0.0, 1000.0)
prob_a = ODEProblem(sys_a, merge(u0_a, params), tspan)
prob_b = ODEProblem(sys_b, merge(u0_b, params), tspan)

sol_a = solve(prob_a, Tsit5(); reltol=1.0e-8, abstol=1.0e-10)
sol_b = solve(prob_b, Tsit5(); reltol=1.0e-8, abstol=1.0e-10)

#

duffing_lab = simplify(trigonometric_form(duffing))
eqs_lab = meanfield([a], duffing_lab, [a]; rates=[γ], order=1, iv=t)
sys_lab = mtkcompile(System(eqs_lab; name=:lab))
prob_lab = ODEProblem(sys_lab, merge(u0_a, params), tspan)

sol_lab = solve(prob_lab, Tsit5(); reltol=1.0e-8, abstol=1.0e-10)

# Compare the slowly varying x envelope; the lab-frame carrier obscures the benchmark.

using Plots: plot, plot!

ω_value = params[ω]
ω₀_value = params[ω₀]
t_plot = sol_lab.t
α_a_rot = get_solution(sol_a, a, eqs_a).(t_plot)
α_lab_rot = exp.(1im .* ω_value .* t_plot) .* get_solution(sol_lab, a, eqs_lab).(t_plot)
β_b_rot = get_solution(sol_b, b, eqs_b).(t_plot)
x_a_envelope = sqrt(2 / ω₀_value) .* abs.(α_a_rot)
x_lab_envelope = sqrt(2 / ω₀_value) .* abs.(α_lab_rot)
x_b_envelope = sqrt(2 / ω_value) .* abs.(β_b_rot)

p_envelope = plot(
  t_plot, x_a_envelope; label="A: effective", xlabel="t", ylabel="x envelope"
)
plot!(p_envelope, t_plot, x_lab_envelope; label="Lab dynamics")
plot!(p_envelope, t_plot, x_b_envelope; label="B: effective")

p_error = plot(
  t_plot,
  abs.(x_a_envelope .- x_lab_envelope);
  label="A error",
  xlabel="t",
  ylabel="absolute envelope error",
  yscale=:log10,
)
plot!(p_error, t_plot, abs.(x_b_envelope .- x_lab_envelope); label="B error")

plot(p_envelope, p_error; layout=(2, 1))

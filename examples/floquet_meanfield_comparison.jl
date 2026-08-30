# # Counting pump photons: a mean-field comparison

using FloquetExpansions

h_a = FockSpace(:cavity)

@qnumbers a::Destroy(h_a)

@variables ω₀::Real ω::Real α::Real F::Real γ::Real
@variables t::Real

x_a = (a + a') / sqrt(2 * ω₀)
p_a = im * sqrt(ω₀ / 2) * (a' - a)
function duffing_hamiltonian(x, p)
  return ω₀ * (p^2 + x^2) // 2 + α * x^4 // 4 - F * x * exponential_form(cos(ω * t))
end

duffing = duffing_hamiltonian(x_a, p_a)

#

rot_a = Rotation(a, ω * t, t)

H_a = simplify(transform(duffing, rot_a))

#

Q = harmonics(H_a, ω, t)

#

rwa_a = floquet_expansion(Q, VanVleck(), 1)
simplify(effective_hamiltonian(rwa_a))

using OrdinaryDiffEq: Tsit5, solve
using QuantumCumulants
using Plots: plot, plot!

hc_a = FockSpace(:a_mode)
hc_b = FockSpace(:b_mode)
@qnumbers a::Destroy(hc_a) b::Destroy(hc_b)

function rotating_duffing(c, wc, omega0, alpha, force, omega, t)
  x = (c + c') / sqrt(2 * wc)
  A = (wc^2 + omega0^2) / (2 * wc)
  B = (omega0^2 - wc^2) / (4 * wc)
  Hstatic = A * (c' * c) + B * (c^2 + c'^2) + alpha * x^4 / 4

  # Transform the static Hamiltonian directly; insert the drive as exact Fourier
  # components so the harmonic collector never sees cos(...) * exp(...).
  hstatic = harmonics(SQA.simplify(transform(Hstatic, Rotation(c, omega * t, t))), omega, t)
  scale = -force / (2 * sqrt(2 * wc))
  components = Dict{Int,SQA.QAdd}(m => q for (m, q) in hstatic.components)
  components[0] = components[0] + scale * (c + c')
  components[2] = get(components, 2, zero(SQA.QAdd)) + scale * c
  components[-2] = get(components, -2, zero(SQA.QAdd)) + scale * c'
  return PeriodicOperator(components, omega)
end

function compile_meanfield(H, c, rate, t, name)
  eqs = meanfield([c], H, [c]; rates=[rate], order=1, iv=t)
  simplify!(eqs)
  sys = mtkcompile(System(eqs; name=name))
  u0 = initial_values(eqs; defaults=Dict(average(c) => 0.0 + 0.0im))
  return (; sys, u0)
end

function solve_meanfield(model, params, t_end)
  prob = ODEProblem(model.sys, merge(model.u0, params), (0.0, t_end))
  return solve(prob, Tsit5(); reltol=1.0e-8, abstol=1.0e-10)
end

function first_harmonic_amplitude(sol, omega_value, omega0_value, t_end)
  period = 2pi / omega_value
  grid = collect(range(t_end - 8period, t_end; length=1601))
  values = [sqrt(2 / omega0_value) * real(sol(time)[1]) for time in grid]
  dt = grid[2] - grid[1]
  duration = grid[end] - grid[1]
  cosine = cos.(omega_value .* grid)
  sine = sin.(omega_value .* grid)
  cosine_integral =
    dt * (sum(values .* cosine) - (values[1] * cosine[1] + values[end] * cosine[end]) / 2)
  sine_integral =
    dt * (sum(values .* sine) - (values[1] * sine[1] + values[end] * sine[end]) / 2)
  return hypot(2cosine_integral / duration, 2sine_integral / duration)
end

omega0_value = 1.0
alpha_value = -0.02
force_value = 0.08
gamma_value = 0.02
t_end = 600.0
orders = 1:2

H_a = rotating_duffing(a, omega0, omega0, alpha, force, omega, t)
H_b = rotating_duffing(b, omega, omega0, alpha, force, omega, t)

expansions_a = [floquet_expansion(H_a, VanVleck(), order) for order in orders]
expansions_b = [floquet_expansion(H_b, VanVleck(), order) for order in orders]

models_a = [
  compile_meanfield(
    effective_hamiltonian(expansion), a, gamma, t, Symbol(:a_order_, order)
  ) for (order, expansion) in zip(orders, expansions_a)
]
models_b = [
  compile_meanfield(
    effective_hamiltonian(expansion), b, gamma, t, Symbol(:b_order_, order)
  ) for (order, expansion) in zip(orders, expansions_b)
]

x_lab = (a + a') / sqrt(2 * omega0)
H_lab = omega0 * (a' * a) + alpha * x_lab^4 / 4 - force * x_lab * cos(omega * t)
model_full = compile_meanfield(H_lab, a, gamma, t, :full_duffing)

function response_at(detuning, model_full, models_a, models_b)
  omega_value = omega0_value + detuning
  params = Dict(
    omega0 => omega0_value,
    omega => omega_value,
    alpha => alpha_value,
    force => force_value,
    gamma => gamma_value,
  )
  full = solve_meanfield(model_full, params, t_end)
  solutions_a = [solve_meanfield(model, params, t_end) for model in models_a]
  solutions_b = [solve_meanfield(model, params, t_end) for model in models_b]
  reference = first_harmonic_amplitude(full, omega_value, omega0_value, t_end)
  amplitudes_a = [
    sqrt(2 / omega0_value) * abs(solution.u[end][1]) for solution in solutions_a
  ]
  amplitudes_b = [
    sqrt(2 / omega_value) * abs(solution.u[end][1]) for solution in solutions_b
  ]
  return (; reference, amplitudes_a, amplitudes_b)
end

detunings = collect(range(0.05, 0.30; length=6))
responses = [
  response_at(detuning, model_full, models_a, models_b) for detuning in detunings
]
reference = [response.reference for response in responses]
amplitudes_a = [
  [response.amplitudes_a[order] for response in responses] for order in eachindex(orders)
]
amplitudes_b = [
  [response.amplitudes_b[order] for response in responses] for order in eachindex(orders)
]

palette = [:darkorange, :royalblue]
response_plot = plot(
  detunings,
  reference;
  label="full mean field",
  color=:black,
  linewidth=3,
  marker=:circle,
  markersize=4,
  xlabel="detuning omega - omega0",
  ylabel="fundamental amplitude",
  title="Driven Duffing response",
  legend=:outerright,
  grid=true,
  framestyle=:box,
)
for (index, order) in enumerate(orders)
  plot!(
    response_plot,
    detunings,
    amplitudes_a[index];
    label="a basis, order $order",
    color=palette[index],
    linestyle=:dash,
    linewidth=2,
  )
  plot!(
    response_plot,
    detunings,
    amplitudes_b[index];
    label="b basis, order $order",
    color=palette[index],
    linewidth=2,
  )
end

error_plot = plot(;
  xlabel="detuning omega - omega0",
  ylabel="absolute amplitude error",
  yscale=:log10,
  legend=:outerright,
  grid=true,
  framestyle=:box,
)
for (index, order) in enumerate(orders)
  plot!(
    error_plot,
    detunings,
    abs.(amplitudes_a[index] .- reference);
    label="a basis, order $order",
    color=palette[index],
    linestyle=:dash,
    linewidth=2,
  )
  plot!(
    error_plot,
    detunings,
    abs.(amplitudes_b[index] .- reference);
    label="b basis, order $order",
    color=palette[index],
    linewidth=2,
  )
end

comparison_plot = plot(
  response_plot,
  error_plot;
  layout=(2, 1),
  link=:x,
  size=(1100, 800),
  plot_title="Counting pump photons improves the mean-field Floquet expansion",
)

rms_error_a = [
  sqrt(sum((amplitudes_a[index] .- reference) .^ 2) / length(reference)) for
  index in eachindex(orders)
]
rms_error_b = [
  sqrt(sum((amplitudes_b[index] .- reference) .^ 2) / length(reference)) for
  index in eachindex(orders)
]
println("RMS errors for a basis: ", rms_error_a)
println("RMS errors for b basis: ", rms_error_b)

comparison_plot

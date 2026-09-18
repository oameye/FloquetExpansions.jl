using BenchmarkTools: @benchmarkable
using Symbolics: @variables

const FE_BVV_INTERNAL = FloquetExpansions

function bvvi_benchmark_workloads()
  pauli = PauliSpace(:bvvi_benchmark_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω_bvvi_q::Real t_bvvi_q::Real Ω_bvvi_q::Real Δ_bvvi_q::Real
  qubit = Δ_bvvi_q * σz + Ω_bvvi_q * cos(ω_bvvi_q * t_bvvi_q) * σx

  fock = FockSpace(:bvvi_benchmark_kerr)
  a = Destroy(fock, :a)
  @variables ω_bvvi_k::Real t_bvvi_k::Real Δ_bvvi_k::Real
  @variables K_bvvi_k::Real ε_bvvi_k::Real
  number = a' * a
  kerr =
    Δ_bvvi_k * number +
    K_bvvi_k * a'^2 * a^2 +
    ε_bvvi_k * cos(ω_bvvi_k * t_bvvi_k) * (a + a')

  return (
    "Driven qubit" => (qubit, ω_bvvi_q, t_bvvi_q),
    "Driven Kerr resonator" => (kerr, ω_bvvi_k, t_bvvi_k),
  )
end

bvvi_bench_product(left, right) = left * right
bvvi_bench_inverse(harmonic) = 1 // harmonic
bvvi_bench_simplify(value) = FE_BVV_INTERNAL.SQA.simplify(value)

function bvvi_benchmark_context(H, ω, t, order)
  generator = FE_BVV_INTERNAL.harmonics(FE_BVV_INTERNAL.qadd(H), ω, t)
  components = getfield(generator, :components)
  plan = FE_BVV_INTERNAL.compile_bloch_projection_plan(collect(keys(generator)), order)
  bloch = FE_BVV_INTERNAL.evaluate_bloch_projection_plan(
    plan,
    components;
    product=bvvi_bench_product,
    inverse_weight=bvvi_bench_inverse,
    zero_component=getfield(generator, :zero_component),
    simplifier=bvvi_bench_simplify,
  )
  return generator, plan, bloch
end

function bvvi_benchmark_reconstruction(generator, plan, bloch)
  return FE_BVV_INTERNAL.bloch_van_vleck_reconstruction(
    plan,
    bloch;
    product=bvvi_bench_product,
    zero_component=getfield(generator, :zero_component),
    simplifier=bvvi_bench_simplify,
  )
end

function benchmark_bloch_van_vleck_internal!(suite)
  for (label, (H, ω, t)) in bvvi_benchmark_workloads(), order in 2:4
    generator, plan, bloch = bvvi_benchmark_context(H, ω, t, order)
    suite["Internal Bloch Van Vleck"][label]["order $order"]["reconstruction"] =
      @benchmarkable bvvi_benchmark_reconstruction($generator, $plan, $bloch)
  end
  return nothing
end

using Test
using FloquetExpansions
using Symbolics: @variables

include(joinpath(@__DIR__, "helpers", "lyndon_log_evaluator.jl"))

const FE_LLE = FloquetExpansions

lle_product(left, right) = left * right
lle_simplify(value) = FE_LLE.SQA.simplify(value)

function lle_qubit_workload()
  pauli = PauliSpace(:lyndon_log_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω_lle_q::Real t_lle_q::Real Ω_lle_q::Real Δ_lle_q::Real
  H = Δ_lle_q * σz + Ω_lle_q * cos(ω_lle_q * t_lle_q) * σx
  return H, ω_lle_q, t_lle_q
end

function lle_kerr_workload()
  fock = FockSpace(:lyndon_log_kerr)
  a = Destroy(fock, :a)
  @variables ω_lle_k::Real t_lle_k::Real Δ_lle_k::Real K_lle_k::Real ε_lle_k::Real
  number = a' * a
  H =
    Δ_lle_k * number +
    K_lle_k * a'^2 * a^2 +
    ε_lle_k * cos(ω_lle_k * t_lle_k) * (a + a')
  return H, ω_lle_k, t_lle_k
end

function lle_context(H, ω, t, order)
  generator = FE_LLE.harmonics(FE_LLE.qadd(H), ω, t)
  components = Dict(harmonic => generator[harmonic] for harmonic in keys(generator))
  zero_component = generator.zero_component
  bloch_plan = FE_LLE.compile_bloch_projection_plan(keys(components), order)
  bloch = FE_LLE.evaluate_bloch_projection_plan(
    bloch_plan,
    components;
    product=lle_product,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component,
    simplifier=lle_simplify,
  )
  direct = FE_LLE.bloch_van_vleck_reconstruction(
    bloch_plan,
    bloch;
    product=lle_product,
    zero_component,
    simplifier=lle_simplify,
  )
  lyndon = compile_lyndon_log_evaluation_plan(keys(components), order)
  evaluated = evaluate_lyndon_log_plan(
    lyndon,
    components;
    product=lle_product,
    zero_component,
    simplifier=lle_simplify,
  )
  return components, zero_component, direct, lyndon, evaluated
end

function lle_embedding_equal(left, right, zero_component)
  harmonics = union(keys(left), keys(right))
  return all(
    iszero(lle_simplify(get(left, harmonic, zero_component) - get(right, harmonic, zero_component))) for
    harmonic in harmonics
  )
end

@testset "compiled Lyndon DAG reproduces SQA canonical logarithm" begin
  for workload in (lle_qubit_workload(), lle_kerr_workload())
    H, ω, t = workload
    for order in 2:4
      _, zero_component, direct, lyndon, evaluated = lle_context(H, ω, t, order)
      @test length(evaluated) == length(direct.log_embedding)
      @test all(
        lle_embedding_equal(evaluated[n], direct.log_embedding[n], zero_component) for
        n in eachindex(evaluated)
      )
      @test lyndon_backend_products(lyndon) == 2 * length(lyndon.brackets)
    end
  end
end

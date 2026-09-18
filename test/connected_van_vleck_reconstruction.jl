using Test
using FloquetExpansions
using Symbolics: @variables

include(joinpath(@__DIR__, "helpers", "connected_van_vleck_reconstruction.jl"))

const FE_CVVT = FloquetExpansions

cvvt_product(left, right) = left * right
cvvt_simplify(value) = FE_CVVT.SQA.simplify(value)

function cvvt_embedding_equal(left, right, zero_component; simplifier=identity)
  harmonics = union(keys(left), keys(right))
  return all(
    iszero(
      simplifier(get(left, harmonic, zero_component) - get(right, harmonic, zero_component))
    ) for harmonic in harmonics
  )
end

@testset "connected reconstruction matches dense noncommuting order six" begin
  R = Rational{Int}
  support = [-1, 0, 1]
  components = Dict(-1 => R[0 1; 2 -1], 0 => R[1 2; -1 0], 1 => R[2 -1; 1 1])
  zero_component = zeros(R, 2, 2)
  order = 6
  projection_plan = FE_CVVT.compile_bloch_projection_plan(support, order)
  bloch = FE_CVVT.evaluate_bloch_projection_plan(
    projection_plan,
    components;
    product=(*),
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component,
  )
  direct = FE_CVVT.bloch_van_vleck_reconstruction(
    projection_plan, bloch; product=(*), zero_component
  )
  log_plan = compile_lyndon_log_evaluation_plan(support, order)
  connected = connected_van_vleck_reconstruction(
    projection_plan, bloch, components, log_plan; product=(*), zero_component
  )

  @test connected.static_factor == direct.static_factor
  @test connected.inverse_static_factor == direct.inverse_static_factor
  @test connected.effective == direct.effective
  @test all(
    cvvt_embedding_equal(
      connected.log_embedding[n], direct.log_embedding[n], zero_component
    ) for n in eachindex(connected.log_embedding)
  )
end

function cvvt_qubit_workload()
  pauli = PauliSpace(:connected_vv_qubit)
  σx = Pauli(pauli, :sigma, 1)
  σz = Pauli(pauli, :sigma, 3)
  @variables ω_cvvt_q::Real t_cvvt_q::Real Ω_cvvt_q::Real Δ_cvvt_q::Real
  H = Δ_cvvt_q * σz + Ω_cvvt_q * cos(ω_cvvt_q * t_cvvt_q) * σx
  return H, ω_cvvt_q, t_cvvt_q
end

function cvvt_kerr_workload()
  fock = FockSpace(:connected_vv_kerr)
  a = Destroy(fock, :a)
  @variables ω_cvvt_k::Real t_cvvt_k::Real Δ_cvvt_k::Real K_cvvt_k::Real ε_cvvt_k::Real
  number = a' * a
  H =
    Δ_cvvt_k * number +
    K_cvvt_k * a'^2 * a^2 +
    ε_cvvt_k * cos(ω_cvvt_k * t_cvvt_k) * (a + a')
  return H, ω_cvvt_k, t_cvvt_k
end

@testset "connected reconstruction matches SQA canonical Van Vleck" begin
  for (H, ω, t) in (cvvt_qubit_workload(), cvvt_kerr_workload())
    generator = FE_CVVT.harmonics(FE_CVVT.qadd(H), ω, t)
    components = Dict(harmonic => generator[harmonic] for harmonic in keys(generator))
    zero_component = generator.zero_component

    for order in 2:4
      projection_plan = FE_CVVT.compile_bloch_projection_plan(keys(components), order)
      bloch = FE_CVVT.evaluate_bloch_projection_plan(
        projection_plan,
        components;
        product=cvvt_product,
        inverse_weight=harmonic -> 1 // harmonic,
        zero_component,
        simplifier=cvvt_simplify,
      )
      direct = FE_CVVT.bloch_van_vleck_reconstruction(
        projection_plan,
        bloch;
        product=cvvt_product,
        zero_component,
        simplifier=cvvt_simplify,
      )
      log_plan = compile_lyndon_log_evaluation_plan(keys(components), order)
      connected = connected_van_vleck_reconstruction(
        projection_plan,
        bloch,
        components,
        log_plan;
        product=cvvt_product,
        zero_component,
        simplifier=cvvt_simplify,
      )

      @test length(connected.static_factor) == length(direct.static_factor)
      @test all(
        iszero(cvvt_simplify(connected.static_factor[n] - direct.static_factor[n])) for
        n in eachindex(connected.static_factor)
      )
      @test all(
        iszero(
          cvvt_simplify(
            connected.inverse_static_factor[n] - direct.inverse_static_factor[n]
          ),
        ) for n in eachindex(connected.inverse_static_factor)
      )
      @test all(
        cvvt_embedding_equal(
          connected.log_embedding[n],
          direct.log_embedding[n],
          zero_component;
          simplifier=cvvt_simplify,
        ) for n in eachindex(connected.log_embedding)
      )
      @test all(
        iszero(cvvt_simplify(connected.effective[n] - direct.effective[n])) for
        n in eachindex(connected.effective)
      )
    end
  end
end

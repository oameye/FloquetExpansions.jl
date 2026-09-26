using Test
using FloquetExpansions
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA_CP_POLICY = SecondQuantizedAlgebra

function cp_policy_vanishes(value)
  return iszero(SQA_CP_POLICY.simplify(value))
end

function cp_policy_vanishes(generator::PeriodicGenerator)
  return iszero(SQA_CP_POLICY.simplify(generator))
end

@testset "CP selector policy is static and enabled by default" begin
  hd_raw = @inferred(HoriDeprit(; complete_positive=Val(false)))
  bf_raw = @inferred(BlochFeshbach(; complete_positive=Val(false)))
  hd_cp = @inferred(HoriDeprit(; complete_positive=Val(true)))
  bf_cp = @inferred(BlochFeshbach(; complete_positive=Val(true)))
  hd_default = @inferred(HoriDeprit())
  bf_default = @inferred(BlochFeshbach())

  @test hd_raw isa HoriDeprit
  @test bf_raw isa BlochFeshbach
  @test typeof(hd_default) === typeof(hd_cp)
  @test typeof(bf_default) === typeof(bf_cp)
  @test typeof(hd_default) !== typeof(hd_raw)
  @test typeof(bf_default) !== typeof(bf_raw)
end

@testset "explicit non-CP selectors preserve raw HD/BF equivalence" begin
  space = PauliSpace(:cp_algorithm_policy)
  σx = Pauli(space, :σ, 1)
  σy = Pauli(space, :σ, 2)
  σz = Pauli(space, :σ, 3)
  @variables ω_cp_policy::Real

  H1 = σx + im * σy
  H = PeriodicGenerator(Dict(0 => 1 * σz, 1 => H1, -1 => H1'), ω_cp_policy)
  L = PeriodicGenerator(
    Dict(
      0 => hamiltonian_action(σz) + dissipator(σx),
      1 => hamiltonian_action(σx + σy),
      -1 => dissipator(σy + σz),
    ),
    ω_cp_policy,
  )

  for generator in (H, L)
    hd = floquet_expansion(
      generator, VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))), 3
    )
    bf = floquet_expansion(
      generator, VanVleck(; algorithm=BlochFeshbach(; complete_positive=Val(false))), 3
    )

    for n in 0:2
      @test cp_policy_vanishes(effective_component(bf, n) - effective_component(hd, n))
    end
    for n in 1:2
      @test cp_policy_vanishes(micromotion(bf, n) - micromotion(hd, n))
    end
  end
end

@testset "CP selector leaves Hamiltonian HFE unchanged" begin
  space = PauliSpace(:cp_algorithm_policy_hamiltonian)
  σx = Pauli(space, :σ, 1)
  σz = Pauli(space, :σ, 3)
  @variables ω_cp_policy_h::Real

  H = PeriodicGenerator(Dict(0 => 1 * σz, 1 => 1 * σx, -1 => 1 * σx), ω_cp_policy_h)

  for (algorithm, raw_algorithm) in (
    (HoriDeprit(), HoriDeprit(; complete_positive=Val(false))),
    (BlochFeshbach(), BlochFeshbach(; complete_positive=Val(false))),
  )
    configured = floquet_expansion(H, VanVleck(; algorithm=algorithm), 3)
    raw = floquet_expansion(H, VanVleck(; algorithm=raw_algorithm), 3)
    @test cp_policy_vanishes(effective_generator(configured) - effective_generator(raw))
    @test cp_policy_vanishes(micromotion(configured) - micromotion(raw))
  end
end

@testset "default Liouvillian selectors apply Gram completion" begin
  space = PauliSpace(:cp_algorithm_policy_liouvillian)
  σx = Pauli(space, :σ, 1)
  σz = Pauli(space, :σ, 3)
  @variables ω_cp_policy_l::Real

  L = PeriodicGenerator(Dict(0 => hamiltonian_action(σz) + dissipator(σx)), ω_cp_policy_l)

  for (algorithm, explicit_cp_algorithm, raw_algorithm) in (
    (
      HoriDeprit(),
      HoriDeprit(; complete_positive=Val(true)),
      HoriDeprit(; complete_positive=Val(false)),
    ),
    (
      BlochFeshbach(),
      BlochFeshbach(; complete_positive=Val(true)),
      BlochFeshbach(; complete_positive=Val(false)),
    ),
  )
    configured = floquet_expansion(L, VanVleck(; algorithm=algorithm), 2)
    explicit_cp = floquet_expansion(L, VanVleck(; algorithm=explicit_cp_algorithm), 2)
    raw = floquet_expansion(L, VanVleck(; algorithm=raw_algorithm), 2)
    expected = positive_completion(raw, Gram())

    @test factorization(configured) isa FloquetExpansions.GramFactorization
    @test cp_policy_vanishes(
      effective_generator(configured) - effective_generator(expected)
    )
    @test cp_policy_vanishes(
      effective_generator(explicit_cp) - effective_generator(configured)
    )
    @test cp_policy_vanishes(
      liouvillian(hamiltonian(configured); channels=channels(configured)) -
      effective_generator(configured),
    )
    for n in 0:1
      @test cp_policy_vanishes(
        effective_component(configured, n) - effective_component(raw, n)
      )
    end
    @test cp_policy_vanishes(micromotion(configured) - micromotion(raw))
  end
end

@testset "CP selector accepts only static Boolean values" begin
  @test_throws ArgumentError HoriDeprit(; complete_positive=Val(:invalid))
  @test_throws ArgumentError BlochFeshbach(; complete_positive=Val(:invalid))
end

using Test
using FloquetExpansions
using LinearAlgebra: norm
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables

const SQA = SecondQuantizedAlgebra

space = NLevelSpace(:native_gksl_reference, 3)
σ11 = Transition(space, :σ, 1, 1)
σ22 = Transition(space, :σ, 2, 2)
σ33 = Transition(space, :σ, 3, 3)
σ12 = Transition(space, :σ, 1, 2)
σ23 = Transition(space, :σ, 2, 3)
σ31 = Transition(space, :σ, 3, 1)
σ13 = Transition(space, :σ, 1, 3)
σ21 = Transition(space, :σ, 2, 1)
σ32 = Transition(space, :σ, 3, 2)
@variables ω_ref::Real

H0_ref = 2 * σ11 - σ22 + (3 // 2) * σ33 + σ12 + σ12'
H1_ref = σ12 + 2 * σ23 + im * σ31
H2_ref = 2 * σ13 - σ21 + im * σ32
H_ref = PeriodicGenerator(
  Dict(0 => H0_ref, 1 => H1_ref, -1 => H1_ref', 2 => H2_ref, -2 => H2_ref'), ω_ref
)
H_map_ref = PeriodicGenerator(
  Dict(harmonic => hamiltonian_action(H_ref[harmonic]) for harmonic in keys(H_ref)), ω_ref
)

@testset "native GKSL reference: Hamiltonian-sector map micromotion through second order" begin
  expansion = floquet_expansion(H_map_ref, VanVleck(), 3)
  micromotion2 = micromotion(expansion, 2)
  nonzero_harmonics = filter(!=(0), collect(keys(H_map_ref)))
  generated_bound = 2 * maximum(abs, nonzero_harmonics)
  @test issubset(keys(micromotion2), (-generated_bound):generated_bound)

  for harmonic in (-generated_bound):generated_bound
    iszero(harmonic) && continue
    expected = (1 // harmonic^2) * SQA.commutator(H_map_ref[harmonic], H_map_ref[0])
    for inner_harmonic in nonzero_harmonics
      inner_harmonic == harmonic && continue
      coefficient = 1 // (2 * harmonic * inner_harmonic)
      expected +=
        coefficient *
        SQA.commutator(H_map_ref[inner_harmonic], H_map_ref[harmonic - inner_harmonic])
    end
    @test iszero(SQA.simplify(ω_ref^2 * micromotion2[harmonic] - expected))
  end
end

R_ref = PeriodicGenerator(
  Dict(
    0 => dissipator(σ23 + im * σ31),
    1 => dissipator(σ12 + σ23),
    -1 => dissipator(σ31 + σ12),
    2 => dissipator(σ13 - im * σ23),
    -2 => dissipator(σ21 + σ32),
    3 => dissipator(σ12 + im * σ31 + σ23),
    -3 => dissipator(σ21 - im * σ13 + σ32),
  ),
  ω_ref,
)

@testset "native GKSL reference: ordered BCH static correction sign" begin
  common_harmonics = filter(!=(0), collect(keys(H_map_ref)))
  cross_average = zero(H_map_ref[0])
  static_correction = zero(H_map_ref[0])

  for harmonic in common_harmonics
    G_harmonic = (im // harmonic) * H_map_ref[harmonic]
    F_minus_harmonic = (im // harmonic) * R_ref[-harmonic]
    cross_average += SQA.commutator(G_harmonic, F_minus_harmonic)
    static_correction +=
      (1 // (2 * harmonic^2)) * SQA.commutator(H_map_ref[harmonic], R_ref[-harmonic])
  end

  @test iszero(SQA.simplify((-1 // 2) * cross_average - static_correction))
  @test !iszero(SQA.simplify((1 // 2) * cross_average - static_correction))
end

@testset "native GKSL reference: one-dissipator second-order sector identities" begin
  h_harmonics = filter(!=(0), collect(keys(H_map_ref)))
  r_harmonics = filter(!=(0), collect(keys(R_ref)))
  three_harmonic_bound = max(2 * maximum(abs, h_harmonics), maximum(abs, r_harmonics))
  three_harmonics = filter(!=(0), collect((-three_harmonic_bound):three_harmonic_bound))
  H0 = H_map_ref[0]
  R0 = R_ref[0]

  C_H0 = zero(H0)
  C_R0 = zero(H0)
  B_R = zero(H0)
  V_R_a = zero(H0)

  for m in h_harmonics
    C_H0 -= (1 // m^2) * SQA.commutator(SQA.commutator(H_map_ref[m], H0), R_ref[-m])
    C_R0 +=
      (1 // (2 * m^2)) * SQA.commutator(H_map_ref[m], SQA.commutator(H_map_ref[-m], R0))
    B_R += (1 // (2 * m^2)) * SQA.commutator(H_map_ref[m], R_ref[-m])

    V_R_a -=
      (1 // (2 * m^2)) * (
        SQA.commutator(R_ref[-m], SQA.commutator(H0, H_map_ref[m])) +
        SQA.commutator(H_map_ref[-m], SQA.commutator(R0, H_map_ref[m])) +
        SQA.commutator(H_map_ref[-m], SQA.commutator(H0, R_ref[m]))
      )
  end

  @test iszero(SQA.simplify(V_R_a - C_H0 - C_R0 - SQA.commutator(B_R, H0)))

  C_3h = zero(H0)
  V_R_b = zero(H0)

  # The dummy indices of the three-harmonic identities can be generated sums outside the
  # original Hamiltonian support. Sum over a closed finite domain instead of the input support.
  for m in three_harmonics, n in three_harmonics
    if n != m
      C_3h -=
        (1 // (2 * m * n)) *
        SQA.commutator(SQA.commutator(H_map_ref[n], H_map_ref[m - n]), R_ref[-m])

      V_R_b -=
        (1 // (3 * n * m)) * (
          SQA.commutator(R_ref[-n], SQA.commutator(H_map_ref[n - m], H_map_ref[m])) +
          SQA.commutator(H_map_ref[-n], SQA.commutator(R_ref[n - m], H_map_ref[m])) +
          SQA.commutator(H_map_ref[-n], SQA.commutator(H_map_ref[n - m], R_ref[m]))
        )
    end

    if m + n != 0
      C_3h -=
        (1 // (2 * m * n)) *
        SQA.commutator(H_map_ref[m], SQA.commutator(H_map_ref[n], R_ref[-(m + n)]))
    end
  end

  @test iszero(SQA.simplify(V_R_b - C_3h))

  V_R = V_R_a + V_R_b
  R_native = C_H0 + C_R0 + C_3h
  @test iszero(SQA.simplify(V_R - R_native - SQA.commutator(B_R, H0)))

  # The order-2 effective component of H + εR is a cubic polynomial in ε, so this four-point
  # stencil returns its linear coefficient exactly: the one-dissipator van Vleck sector.
  second_order(ε) =
    ω_ref^2 *
    effective_component(floquet_expansion(H_map_ref + ε * R_ref, VanVleck(), 3), 2)
  one_dissipator_van_vleck =
    (1 // 6) *
    (6 * second_order(1) - 2 * second_order(-1) - 3 * second_order(0) - second_order(2))
  @test iszero(SQA.simplify(one_dissipator_van_vleck - V_R))
end

matrix_commutator_ref(A, B) = A * B - B * A

@testset "native GKSL reference: one-dissipator generic matrix sanity oracle" begin
  zero_matrix = zeros(ComplexF64, 2, 2)
  H_matrix = Dict(
    0 => ComplexF64[0.3 0.7 - 0.2im; -0.4 + 0.1im -0.6],
    1 => ComplexF64[0.2im 0.8; -0.3 0.4 - 0.1im],
    -1 => ComplexF64[0.5 -0.2im; 0.6im -0.1],
    2 => ComplexF64[-0.4im 0.3; 0.7 + 0.2im 0.5],
    -2 => ComplexF64[0.1 0.9im; -0.6 0.2 - 0.3im],
  )
  R_matrix = Dict(
    0 => ComplexF64[-0.2 0.4im; 0.5 0.3],
    1 => ComplexF64[0.6 0.1 - 0.2im; -0.4im -0.5],
    -1 => ComplexF64[0.3im -0.7; 0.2 + 0.1im 0.4],
    2 => ComplexF64[-0.1 0.5; 0.8im 0.2],
    -2 => ComplexF64[0.7 -0.3im; -0.2 0.6im],
    3 => ComplexF64[0.4 + 0.2im -0.5; 0.3im -0.6],
    -3 => ComplexF64[-0.2im 0.8; -0.1 + 0.4im 0.7],
  )
  harmonic(H, m) = get(H, m, zero_matrix)
  h_harmonics = [-2, -1, 1, 2]
  r_harmonics = [-3, -2, -1, 1, 2, 3]
  three_harmonic_bound = max(2 * maximum(abs, h_harmonics), maximum(abs, r_harmonics))
  three_harmonics = filter(!=(0), collect((-three_harmonic_bound):three_harmonic_bound))
  H0 = H_matrix[0]
  R0 = R_matrix[0]

  C_H0 = copy(zero_matrix)
  C_R0 = copy(zero_matrix)
  B_R = copy(zero_matrix)
  V_R_a = copy(zero_matrix)

  for m in h_harmonics
    Hm = harmonic(H_matrix, m)
    Hminus = harmonic(H_matrix, -m)
    Rm = harmonic(R_matrix, m)
    Rminus = harmonic(R_matrix, -m)

    C_H0 -= matrix_commutator_ref(matrix_commutator_ref(Hm, H0), Rminus) / m^2
    C_R0 += matrix_commutator_ref(Hm, matrix_commutator_ref(Hminus, R0)) / (2 * m^2)
    B_R += matrix_commutator_ref(Hm, Rminus) / (2 * m^2)
    V_R_a -=
      (
        matrix_commutator_ref(Rminus, matrix_commutator_ref(H0, Hm)) +
        matrix_commutator_ref(Hminus, matrix_commutator_ref(R0, Hm)) +
        matrix_commutator_ref(Hminus, matrix_commutator_ref(H0, Rm))
      ) / (2 * m^2)
  end

  C_3h = copy(zero_matrix)
  V_R_b = copy(zero_matrix)

  for m in three_harmonics, n in three_harmonics
    if n != m
      C_3h -=
        matrix_commutator_ref(
          matrix_commutator_ref(harmonic(H_matrix, n), harmonic(H_matrix, m - n)),
          harmonic(R_matrix, -m),
        ) / (2 * m * n)

      V_R_b -=
        (
          matrix_commutator_ref(
            harmonic(R_matrix, -n),
            matrix_commutator_ref(harmonic(H_matrix, n - m), harmonic(H_matrix, m)),
          ) +
          matrix_commutator_ref(
            harmonic(H_matrix, -n),
            matrix_commutator_ref(harmonic(R_matrix, n - m), harmonic(H_matrix, m)),
          ) +
          matrix_commutator_ref(
            harmonic(H_matrix, -n),
            matrix_commutator_ref(harmonic(H_matrix, n - m), harmonic(R_matrix, m)),
          )
        ) / (3 * n * m)
    end

    if m + n != 0
      C_3h -=
        matrix_commutator_ref(
          harmonic(H_matrix, m),
          matrix_commutator_ref(harmonic(H_matrix, n), harmonic(R_matrix, -(m + n))),
        ) / (2 * m * n)
    end
  end

  R_native = C_H0 + C_R0 + C_3h
  V_R = V_R_a + V_R_b
  residual = V_R - R_native - matrix_commutator_ref(B_R, H0)

  @test norm(V_R_a - C_H0 - C_R0 - matrix_commutator_ref(B_R, H0)) <= 5.0e-13
  @test norm(V_R_b - C_3h) <= 5.0e-13
  @test norm(residual) <= 5.0e-13
end

@testset "native GKSL reference: RR sector begins at quadratic dissipative strength" begin
  R_osc = PeriodicGenerator(
    Dict(harmonic => R_ref[harmonic] for harmonic in keys(R_ref) if !iszero(harmonic)),
    ω_ref,
  )
  expansion = floquet_expansion(R_osc, VanVleck(), 2)
  expected = zero(R_ref[0])

  for harmonic in keys(R_osc)
    expected += (im // harmonic) * compose(R_osc[-harmonic], R_osc[harmonic])
  end

  @test iszero(SQA.simplify(ω_ref * effective_component(expansion, 1) - expected))
  @test !iszero(SQA.simplify(expected))
end

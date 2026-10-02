using Test
using LinearAlgebra: Diagonal, norm
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions
const Q = Complex{Rational{BigInt}}
const R = Rational{BigInt}

frame_coordinates(rep, X) = (rep.inverse_frame * vec(X))[2:end]

# Driven qubit with static loss written in a possibly rescaled basis |k) = sqrt(g_k)|k>, where the
# operator matrices are D^{-1} X D with D = diag(1, r) and the Hilbert metric is g = (1, r^2).
function exact_qubit(r)
  rep = FE.ExactLiouvilleRepresentation{Q}(R[1, r ^ 2])
  D = Diagonal(Q[1, r])
  rescale(X) = inv(D) * X * D
  σx = rescale(Q[0 1; 1 0])
  σz = rescale(Q[1 0; 0 -1])
  σm = rescale(Q[0 0; 1 0])
  c = frame_coordinates(rep, σm)
  zero3 = zeros(Q, 3, 3)
  L = Dict{Int,Matrix{Q}}(
    0 => FE.native_gksl(rep, (3 // 10) * σz, (3 // 5) * c * c'),
    1 => FE.native_gksl(rep, (2 // 5) * σx, zero3),
    -1 => FE.native_gksl(rep, (2 // 5) * σx, zero3),
  )
  return rep, L, reshape(c, :, 1), kron(Matrix(D), Matrix(inv(D)))
end

@testset "exact representation reconstructs a GKSL generator in a rescaled basis" begin
  g = R[1, 1, 2]
  rep = FE.ExactLiouvilleRepresentation{Q}(g)
  a = Q[0 1 0; 0 0 2; 0 0 0]
  @test FE.hilbert_adjoint(g, a) == Q[0 0 0; 1 0 0; 0 1 0]
  H = Q[1//3 1//2 0; 1//2 -1 2im; 0 -1im 1//5]
  H = (H + FE.hilbert_adjoint(g, H)) / 2
  c = frame_coordinates(rep, a)
  C = (3 // 4) * c * c'
  L = FE.native_gksl(rep, H, C)
  @test FE.native_kossakowski(rep, L) == C
  Hread = FE.native_hamiltonian(rep, L)
  traceless(X) = X - (sum(X[k, k] for k in 1:3) // 3) * Matrix{Q}(FE.LinearAlgebra.I, 3, 3)
  @test Hread == traceless(H)
  @test FE.native_gksl(rep, Hread, C) == L
  # Every gauge generator is trace annihilating: tr S(ρ) = 0 for every matrix unit ρ.
  tracevec = vec(Matrix{Q}(FE.LinearAlgebra.I, 3, 3))
  @test all(iszero(transpose(tracevec) * G) for G in rep.gauge)
  @test length(rep.gauge) == 3^4 - 3^2
end

@testset "exact native recurrence: BF equals HD exactly and matches the float core" begin
  rep, L, leading, _ = exact_qubit(1)
  for inverse in
      (FE.NoHomologicalInverse(), FE.ChargeGradedInverse(rep, L[0], leading, 1e-9))
    bf = FE.native_recurrence(
      FE.BlochFeshbach(), rep, inverse, L, leading, Q[3 // 5], 4, 1e-8
    )
    hd = FE.native_recurrence(FE.HoriDeprit(), rep, inverse, L, leading, Q[3 // 5], 4, 1e-8)
    @test bf isa FE.NativeRecurrence{Matrix{Q},Matrix{Q},Q}
    @test bf.E == hd.E
    @test bf.H == hd.H
    for order in 0:4
      @test FE.native_kossakowski(rep, bf.E[order + 1]) ==
        FE.gram_coefficient(bf.channels, order, 3)
    end
  end

  float = FE.DenseLiouvilleRepresentation(2)
  σx = ComplexF64[0 1; 1 0]
  σz = ComplexF64[1 0; 0 -1]
  σm = ComplexF64[0 0; 1 0]
  c = ComplexF64[FE.LinearAlgebra.tr(F' * σm) for F in float.basis]
  zero3 = zeros(ComplexF64, 3, 3)
  Lf = Dict{Int,Matrix{ComplexF64}}(
    0 => FE.native_gksl(float, 0.3 * σz, 0.6 * c * c'),
    1 => FE.native_gksl(float, 0.4 * σx, zero3),
    -1 => FE.native_gksl(float, 0.4 * σx, zero3),
  )
  leadf = reshape(c, :, 1)
  ref = FE.native_recurrence(
    FE.BlochFeshbach(),
    float,
    FE.NoHomologicalInverse(),
    Lf,
    leadf,
    ComplexF64[0.6],
    4,
    1e-9,
  )
  exact = FE.native_recurrence(
    FE.BlochFeshbach(), rep, FE.NoHomologicalInverse(), L, leading, Q[3 // 5], 4, 1e-8
  )
  @test all(norm(ComplexF64.(exact.E[k]) - ref.E[k]) <= 1e-12 for k in 1:5)
  @test all(norm(ComplexF64.(exact.S[k]) - ref.S[k]) <= 1e-12 for k in 1:4)
end

@testset "exact native recurrence is covariant under a Hilbert metric" begin
  rep1, L1, lead1, _ = exact_qubit(1)
  rep2, L2, lead2, M = exact_qubit(2)
  transform(X) = M * X * inv(M)
  @test transform(L1[0]) == L2[0]
  for inverse in (
    (FE.NoHomologicalInverse(), FE.NoHomologicalInverse()),
    (
      FE.ChargeGradedInverse(rep1, L1[0], lead1, 1e-9),
      FE.ChargeGradedInverse(rep2, L2[0], lead2, 1e-9),
    ),
  )
    r1 = FE.native_recurrence(
      FE.BlochFeshbach(), rep1, inverse[1], L1, lead1, Q[3 // 5], 4, 1e-8
    )
    r2 = FE.native_recurrence(
      FE.BlochFeshbach(), rep2, inverse[2], L2, lead2, Q[3 // 5], 4, 1e-8
    )
    @test all(transform(r1.E[k]) == r2.E[k] for k in 1:5)
    @test all(transform(r1.S[k]) == r2.S[k] for k in 1:4)
  end
end

@testset "exact native slot rejects an indefinite canonical residual" begin
  rep = FE.ExactLiouvilleRepresentation{Q}(2)
  L0 = FE.native_gksl(rep, Q[1 0; 0 -1], zeros(Q, 3, 3))
  residual = FE.native_gksl(rep, zeros(Q, 2, 2), Q[1 0 0; 0 -1 0; 0 0 0])
  @test_throws ArgumentError FE.native_static_step(
    rep, FE.NoHomologicalInverse(), L0, residual, zeros(Q, 3, 3), zeros(Q, 3, 0), Q[], 1e-8
  )
end

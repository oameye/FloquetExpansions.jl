using Test
using LinearAlgebra: Diagonal, I
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions
const Q = Complex{Rational{BigInt}}
const R = Rational{BigInt}

# Four Kossakowski directions: e1, e2 active, e3, e4 dark. The gauge images have dark-dark parts
# orthogonal to the newborn form diag(0, 0, 2, 0), so the canonical least-squares section removes
# them exactly and leaves a PSD residual with one channel born along e3.
function exact_fixture(weight)
  B = Q[1 0; 1//2 1; 0 0; 0 0]
  weights = R[weight, 1]
  hermitian(X) = (X + X') / 2
  g1 = hermitian(Q[1 2 0 1; 0 0 1im 0; 0 0 0 1; 0 0 0 0])
  g2 = hermitian(Q[0 1//3 1 0; 0 1 0 0; 0 0 0 1im; 0 0 0 0])
  g3 = hermitian(Q[0 0 0 0; 0 0 0 2//5; 0 0 0 0; 0 0 0 1])
  images = [g1, g2, g3]
  known = hermitian(Q[1//4 0 0 0; 1im 1//2 0 0; 0 0 0 0; 0 0 0 0])
  Y0 = Q[1//3 0; 0 -1im; 1//5 2; 0 1//7]
  Bw = B * Diagonal(Q.(weights))
  newborn = Q[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 0]
  residual = known + Bw * Y0' + Y0 * Bw' + newborn - (3 * g1 - g2 // 2 + 2 * g3)
  return (; B, weights, images, known, residual)
end

tofloat(X) = ComplexF64.(X)

@testset "exact static slot matches the float core in orthonormal coordinates" begin
  for weight in (1, 4)
    f = exact_fixture(weight)
    G = Matrix{Q}(I, 4, 4)
    Γ = Matrix{R}(I, 3, 3)
    exact = FE.native_exact_static_solve(
      f.residual, f.known, f.B, f.weights, f.images, G, Γ
    )
    @test exact isa FE.ExactStaticSolution{Q,R}
    @test exact.coordinates == R[3, -1 // 2, 2]
    @test exact.newborn_weights == R[2]
    @test exact.dark_residual == Q[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 0]

    amplitudes = tofloat(f.B * Diagonal(Q.(sqrt.(f.weights))))
    float = FE.native_static_solve(
      tofloat(f.residual), tofloat(f.known), amplitudes, tofloat.(f.images), 1e-10
    )
    @test float.coordinates ≈ Float64.(exact.coordinates) atol = 1e-9
    @test float.coefficient ≈ tofloat(exact.coefficient) atol = 1e-9
    born = exact.newborn * Diagonal(Q.(exact.newborn_weights)) * exact.newborn'
    @test float.newborn * float.newborn' ≈ tofloat(born) atol = 1e-9
    # Rate-weighted corrections y_j are the physical amplitude corrections divided by √w_j.
    @test float.correction ≈ tofloat(exact.correction * Diagonal(Q.(sqrt.(f.weights)))) atol =
      1e-9
  end
end

@testset "exact static slot is covariant under a non-orthonormal frame" begin
  f = exact_fixture(1)
  G = Matrix{Q}(I, 4, 4)
  Γ = Matrix{R}(I, 3, 3)
  reference = FE.native_exact_static_solve(
    f.residual, f.known, f.B, f.weights, f.images, G, Γ
  )

  A = Q[1 1//2 0 0; 0 1 1//3 0; 0 0 1 1im; 1//4 0 0 1]
  Ainv = inv(A)
  change(X) = Ainv * X * Ainv'
  transformed = FE.native_exact_static_solve(
    change(f.residual),
    change(f.known),
    Ainv * f.B,
    f.weights,
    change.(f.images),
    A' * G * A,
    Γ,
  )
  @test transformed.coordinates == reference.coordinates
  @test transformed.coefficient == change(reference.coefficient)
  @test transformed.dark_residual == change(reference.dark_residual)
  @test transformed.correction == Ainv * reference.correction
end

@testset "exact static slot is covariant under a gauge reparameterization" begin
  f = exact_fixture(1)
  G = Matrix{Q}(I, 4, 4)
  Γ = Matrix{R}(I, 3, 3)
  reference = FE.native_exact_static_solve(
    f.residual, f.known, f.B, f.weights, f.images, G, Γ
  )

  Rm = R[1 1 0; 0 2 0; 1//2 0 1]
  images = [sum(Rm[i, j] * f.images[i] for i in 1:3) for j in 1:3]
  transformed = FE.native_exact_static_solve(
    f.residual, f.known, f.B, f.weights, images, G, Matrix(transpose(Rm)) * Γ * Rm
  )
  @test transformed.coordinates == Rm \ reference.coordinates
  @test transformed.coefficient == reference.coefficient
end

@testset "exact static slot rejects an indefinite canonical residual" begin
  f = exact_fixture(1)
  G = Matrix{Q}(I, 4, 4)
  Γ = Matrix{R}(I, 3, 3)
  shifted = f.residual - Q[0 0 0 0; 0 0 0 0; 0 0 3 0; 0 0 0 0]
  @test_throws FE.NativePositivityError FE.native_exact_static_solve(
    shifted, f.known, f.B, f.weights, f.images, G, Γ
  )
  @test_throws DimensionMismatch FE.native_exact_static_solve(
    f.residual, f.known, f.B, R[1], f.images, G, Γ
  )
end

@testset "explicit exact kernels agree with the LinearAlgebra operations they replace" begin
  A = Q[1 2im 0; 1//2 -1 3; 0 1im 1//3]
  B = Q[2 0 1im; 1 1//5 0; -1 2 1]
  x = Q[1, 1im, -1 // 2]
  @test FE.exact_mul(A, B) == A * B
  @test FE.exact_mul(A, x) == A * x
  @test FE.exact_mul(A, B, A) == A * B * A
  @test FE.exact_mul(A[:, 1:2], B[1:2, :]) == A[:, 1:2] * B[1:2, :]
  @test FE.exact_adjoint(A) == Matrix(A')
  @test FE.exact_adjoint(A[:, 1:2]) == Matrix(A[:, 1:2]')
  @test FE.exact_diagonal(Q, R[2, 3 // 5]) == Matrix(Diagonal(Q[2, 3 // 5]))
  @test_throws DimensionMismatch FE.exact_mul(A, B[1:2, :])

  # Solves with matrix and vector right-hand sides, including a row interchange.
  P = Q[0 1 1im; 2 0 1; 1 1 1//4]
  @test P * FE.exact_solve(P, B) == B
  @test P * FE.exact_solve(P, x) == x
  @test_throws ArgumentError FE.exact_solve(Q[1 2 3; 2 4 6; 0 1 1], B)
  input = copy(B)
  FE.exact_solve(P, input)
  @test input == B

  # Nullspace and row basis of a rank-two matrix over both the complex and real scalar rings.
  N = Q[1 2 3; 2 4 6; 1im 2im 3im + 1]
  kernel = FE.exact_nullspace(N)
  @test length(kernel) == 1
  @test all(iszero, N * kernel[1])
  @test FE.exact_nullspace(R[1 2 3; 2 4 6]) == [R[-2, 1, 0], R[-3, 0, 1]]
  @test FE.exact_row_basis(N) == [1, 3]
  @test FE.exact_row_basis(R[1 2; 2 4; 0 1; 1 3]) == [1, 3]
end

@testset "exact PSD factor does not modify its input" begin
  P = Q[2 1im; -1im 1]
  copied = copy(P)
  newborn, weights = FE.exact_psd_factor(P, R)
  @test P == copied
  @test newborn * Diagonal(Q.(weights)) * newborn' == P
end

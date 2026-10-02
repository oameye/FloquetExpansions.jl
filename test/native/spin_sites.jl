using Test
using Random: MersenneTwister
using LinearAlgebra: Hermitian, eigmin
using Symbolics: @variables
using FloquetExpansions
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions
const Q = Complex{Rational{BigInt}}
const R = Rational{BigInt}

@variables ω::Real t::Real Δ::Real

# Normal-order a word in the letters :p (S₊), :z (S_z), :m (S₋) with the three PBW rewriting
# rules, independently of the closed-form product under test.
const REWRITES = Dict(
  (:z, :p) => [([:p, :z], 1), ([:p], 1)],
  (:m, :z) => [([:z, :m], 1), ([:m], 1)],
  (:m, :p) => [([:p, :m], 1), ([:z], -2)],
)

function normal_order(letters::Vector{Symbol})
  pending = Dict{Vector{Symbol},Int}(letters => 1)
  done = Dict{NTuple{3,Int},Int}()
  while !isempty(pending)
    current, coefficient = first(pending)
    delete!(pending, current)
    position = findfirst(
      i -> haskey(REWRITES, (current[i], current[i + 1])), 1:(length(current) - 1)
    )
    if position === nothing
      key = (count(==(:p), current), count(==(:z), current), count(==(:m), current))
      done[key] = get(done, key, 0) + coefficient
      continue
    end
    head, tail = current[1:(position - 1)], current[(position + 2):end]
    for (middle, factor) in REWRITES[(current[position], current[position + 1])]
      next = vcat(head, middle, tail)
      pending[next] = get(pending, next, 0) + factor * coefficient
    end
  end
  return Dict(k => v for (k, v) in done if v != 0)
end

letters(m::NTuple{3,Int}) = vcat(fill(:p, m[1]), fill(:z, m[2]), fill(:m, m[3]))

spin_algebra() = FE.OperatorAlgebra{Q}([FE.spin_site(R)])

@testset "PBW product agrees with brute-force reordering" begin
  rng = MersenneTwister(7)
  for _ in 1:150
    left = Tuple(rand(rng, 0:3, 3))
    right = Tuple(rand(rng, 0:3, 3))
    expected = normal_order(vcat(letters(left), letters(right)))
    got = Dict(m => Int(real(c)) for (m, c) in FE.spin_site_product(left, right, Q))
    @test got == expected
  end
end

@testset "PBW product agrees with spin matrices" begin
  function matrices(j)
    dim = Int(2j + 1)
    ms = [j - k for k in 0:(dim - 1)]
    Sz = [i == k ? Float64(ms[i]) : 0.0 for i in 1:dim, k in 1:dim]
    Sp = zeros(dim, dim)
    for k in 2:dim
      Sp[k - 1, k] = sqrt(j * (j + 1) - ms[k] * (ms[k] + 1))
    end
    return Sp, Sz, Matrix(Sp')
  end
  function evaluate(terms, (Sp, Sz, Sm))
    return sum(Float64(real(c)) * Sp^a * Sz^b * Sm^f for ((a, b, f), c) in terms)
  end
  rng = MersenneTwister(11)
  for j in (3 // 2, 2, 5 // 2), _ in 1:25
    S = matrices(j)
    left = Tuple(rand(rng, 0:2, 3))
    right = Tuple(rand(rng, 0:2, 3))
    product = FE.spin_site_product(left, right, Q)
    A = S[1]^left[1] * S[2]^left[2] * S[3]^left[3]
    B = S[1]^right[1] * S[2]^right[2] * S[3]^right[3]
    @test evaluate(product, S) ≈ A * B atol = 1e-9
  end
end

@testset "spin algebra relations, adjoint and degree" begin
  algebra = spin_algebra()
  Sp = FE.algebra_operator(algebra, [([1, 0, 0], 1)])
  Sz = FE.algebra_operator(algebra, [([0, 1, 0], 1)])
  Sm = FE.algebra_operator(algebra, [([0, 0, 1], 1)])
  @test Sp * Sm - Sm * Sp == 2 * Sz
  @test Sz * Sp - Sp * Sz == Sp
  @test Sz * Sm - Sm * Sz == -Sm
  @test adjoint(Sp) == Sm
  @test adjoint(Sz) == Sz
  rng = MersenneTwister(3)
  random_operator() = FE.algebra_operator(
    algebra, [(rand(rng, 0:3, 3), rand(rng, -3:3) + im * rand(rng, -3:3)) for _ in 1:3]
  )
  for _ in 1:20
    X, Y = random_operator(), random_operator()
    @test adjoint(adjoint(X)) == X
    @test adjoint(X * Y) == adjoint(Y) * adjoint(X)
  end
  @test FE.monomial_degree(algebra, [2, 1, 3]) == 6
  @test FE.monomial_charge(algebra, [0], 1, [2, 1, 3]) == [-1]
  @test length(FE.algebra_monomials(algebra, 2)) == 9
end

@testset "mixed boson, spin and level sites keep per-site slots" begin
  algebra = FE.OperatorAlgebra{Q}([
    FE.boson_site(R), FE.spin_site(R), FE.level_site(R[1, 1])
  ])
  @test algebra.offsets == [0, 2, 5]
  @test FE.identity_monomial(algebra) == zeros(Int, 7)
  @test FE.site_slots(algebra, 2) == 3:5
  a = FE.algebra_operator(algebra, [([0, 1, 0, 0, 0, 0, 0], 1)])
  Sp = FE.algebra_operator(algebra, [([0, 0, 1, 0, 0, 0, 0], 1)])
  Sm = FE.algebra_operator(algebra, [([0, 0, 0, 0, 1, 0, 0], 1)])
  Sz = FE.algebra_operator(algebra, [([0, 0, 0, 1, 0, 0, 0], 1)])
  @test a * Sp == Sp * a
  @test Sp * Sm - Sm * Sp == 2 * Sz
end

@testset "SQA lowering and lifting of Spin operators round trip" begin
  h = FockSpace(:c) ⊗ SpinSpace(:S)
  a = Destroy(h, :a)
  Sx, Sy, Sz = (Spin(h, :S, k) for k in 1:3)
  H = Δ * Sz + (2 // 5) * (a + a') * Sx * cos(ω * t)
  G = harmonics(liouvillian(H; channels=(jump(Sx - im * Sy, 4 // 5),)), ω, t)
  lowering = FE.SQALowering{Complex{Rational{Int128}}}(G, Dict(Δ => 1 // 2))
  @test lowering.algebra.offsets == [0, 2]
  lowered(op) = FE.lower_qadd(lowering, FE.qadd(op))
  @test lowered(Sx * Sy - Sy * Sx) == im * lowered(Sz)
  @test lowered(Sy * Sz - Sz * Sy) == im * lowered(Sx)
  @test lowered(Sz * Sx - Sx * Sz) == im * lowered(Sy)
  @test adjoint(lowered(Sx + im * Sy)) == lowered(Sx - im * Sy)
  for op in (Sx, Sy, Sz, Sx * Sy, Sx - im * Sy, a' * Sz * Sx * Sx, (Sx + im * Sy)^3 * Sz)
    Y = lowered(op)
    @test FE.lower_qadd(lowering, FE.lift_operator(lowering, Y)) == Y
  end
  rng = MersenneTwister(5)
  for _ in 1:10
    monomial = [0, 0, rand(rng, 0:3), rand(rng, 0:3), rand(rng, 0:3)]
    Y = FE.algebra_operator(lowering.algebra, [(monomial, 1)])
    @test FE.lower_qadd(lowering, FE.lift_operator(lowering, Y)) == Y
  end
  incomplete = FE.SQALowering{Complex{Rational{Int128}}}(G, Dict{Any,Any}())
  @test_throws ArgumentError FE.lower_generator(incomplete, G)
end

# The Van Vleck reference is computed in Float64, so its lowering is compared with a tolerance.
# Its floats lower to rationals with large denominators that overflow Int128 sums, hence the
# reference is lowered in Rational{BigInt} and both sides are compared entrywise in ComplexF64.
function nearly_same(A::FE.AlgebraSuperoperator, B::FE.AlgebraSuperoperator)
  difference(k) = ComplexF64(get(A.terms, k, 0)) - ComplexF64(get(B.terms, k, 0))
  return all(abs(difference(k)) <= 1e-12 for k in union(keys(A.terms), keys(B.terms)))
end

function driven_spin(detuning)
  h = SpinSpace(:S)
  Sx, Sy, Sz = (Spin(h, :S, k) for k in 1:3)
  Sp = Sx + im * Sy
  Sm = Sx - im * Sy
  H = detuning * Sz + (2 // 5) * (Sp + Sm) * cos(ω * t)
  return H, Sm
end

const SPIN_VALUES = Dict(Δ => 1 // 2)

@testset "native expansion of a driven spin equals Van Vleck" begin
  H, Sm = driven_spin(Δ)
  G = harmonics(liouvillian(H; channels=(jump(Sm, 4 // 5),)), ω, t)
  for N in (2, 3), algorithm in (BlochFeshbach(), HoriDeprit())
    data = FE.native_expansion_data(algorithm, G, N, SPIN_VALUES, 3)
    exact = FE.SQALowering{Q}(G, SPIN_VALUES)
    reference = floquet_expansion(G, VanVleck(; algorithm), N + 1)
    @test all(iszero, data.recurrence.S)
    for n in 0:N
      @test nearly_same(
        data.recurrence.E[n + 1],
        FE.lower_liouvillian(exact, reference.effective_components[n + 1]),
      )
    end
  end
end

@testset "GKSLNormalForm on a driven spin is completely positive" begin
  H, Sm = driven_spin(1 // 2)
  native = floquet_expansion(H, ω, t, GKSLNormalForm(), 3; channels=(jump(Sm, 4 // 5),))
  C = kossakowski(native)
  for value in (3.0, 10.0, 100.0)
    point = Dict(ω => value)
    numeric = [ComplexF64(FE.SQA.to_complex(FE.SQA.substitute(c, point))) for c in C]
    @test eigmin(Hermitian(numeric)) >= -1e-12
  end
end

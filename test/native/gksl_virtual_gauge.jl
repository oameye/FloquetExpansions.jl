using Test
using LinearAlgebra: Hermitian, I, eigmin, eigvals
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using FloquetExpansions
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions
const SQA = SecondQuantizedAlgebra

@variables ω::Real t::Real

function same(A::Liouvillian, B::Liouvillian)
  return FE.liouvillian_iszero(FE.canonical_liouvillian(SQA.simplify(A - B)))
end
same(A, B) = iszero(SQA.simplify(FE.qadd(A) - FE.qadd(B)))

const MODE = FockSpace(:virtual_gauge)
const A = Destroy(MODE, :a)
const LOSS = (jump(A, 4 // 5),)

function kerr(c1=1 // 2, extra=0)
  return c1 * A' * A +
         (3 // 10) * A' * A * A' * A +
         (2 // 5) * (A + A') * cos(ω * t) +
         extra
end

function kerr_expansion(order, algorithm=BlochFeshbach(); H=kerr())
  return floquet_expansion(H, ω, t, GKSLNormalForm(; algorithm), order; channels=LOSS)
end

const BF = kerr_expansion(7)
const HD = kerr_expansion(7, HoriDeprit())

function ladder_matrix(d)
  M = zeros(ComplexF64, d, d)
  for n in 1:(d - 1)
    M[n, n + 1] = sqrt(n)
  end
  return M
end

function dense_operator(q, d, value)
  lowering = ladder_matrix(d)
  M = zeros(ComplexF64, d, d)
  for (term, coeff) in SQA.substitute(q, Dict(ω => value))
    z = SQA.to_num(coeff)
    c = complex(Float64(Symbolics.value(real(z))), Float64(Symbolics.value(imag(z))))
    T = Matrix{ComplexF64}(I, d, d)
    for o in term.ops
      T = T * (isequal(o, A) ? lowering : lowering')
    end
    M .+= c .* T
  end
  return M
end

function hamiltonian_superoperator(H::AbstractMatrix)
  Id = Matrix{ComplexF64}(I, size(H)...)
  return -im * (kron(Id, H) - kron(transpose(H), Id))
end

function dissipator_superoperator(L::AbstractMatrix)
  Id = Matrix{ComplexF64}(I, size(L)...)
  LL = L' * L
  return kron(conj(L), L) - kron(Id, LL) / 2 - kron(transpose(LL), Id) / 2
end

function effective_superoperator(expansion, d, value)
  S = hamiltonian_superoperator(dense_operator(hamiltonian(expansion), d, value))
  for c in channels(expansion)
    rate = real(ComplexF64(SQA.to_complex(SQA.substitute(c.rate, Dict(ω => value)))))
    S += rate * dissipator_superoperator(dense_operator(c.operator, d, value))
  end
  return S
end

function floquet_exponents(d, value, harmonics)
  lowering = ladder_matrix(d)
  number = lowering' * lowering
  static =
    hamiltonian_superoperator(number / 2 + 3 * number * number / 10) +
    4 * dissipator_superoperator(lowering) / 5
  drive = hamiltonian_superoperator(2 * (lowering + lowering') / 5)
  D = d^2
  Q = zeros(ComplexF64, D * (2harmonics + 1), D * (2harmonics + 1))
  for (b, k) in enumerate((-harmonics):harmonics)
    block = ((b - 1) * D + 1):(b * D)
    Q[block, block] = static - im * k * value * I
    if b <= 2harmonics
      Q[block, block .+ D] = drive / 2
      Q[block .+ D, block] = drive / 2
    end
  end
  return eigvals(Q)
end

function slow_distance(exact, approximate)
  return maximum(
    minimum(abs.(exact .- λ)) for λ in approximate if real(λ) > -3 // 2; init=0.0
  )
end

function loglog_slope(values, errors)
  return -(log(errors[end]) - log(errors[1])) / (log(values[end]) - log(values[1]))
end

@testset "driven Kerr reaches the last retained order through a virtual gauge" begin
  @test virtual_orders(BF) == [6]
  @test virtual_orders(HD) == [6]
  @test isempty(virtual_orders(kerr_expansion(5)))
  @test isempty(
    virtual_orders(floquet_expansion(kerr(), ω, t, VanVleck(), 5; channels=LOSS))
  )
  for expansion in (BF, HD)
    @test length(channels(expansion)) == 4
    @test same(
      effective_generator(expansion),
      liouvillian(hamiltonian(expansion); channels=Tuple(channels(expansion))),
    )
    @test !haskey(micromotion(expansion, 6).components, 0)
    @test !isempty(micromotion(expansion, 6).components)
  end
end

@testset "Bloch-Feshbach and Hori-Deprit agree at the virtual order" begin
  for n in 0:6
    @test same(effective_component(BF, n), effective_component(HD, n))
  end
  @test same(hamiltonian(BF), hamiltonian(HD))
  for (c, d) in zip(channels(BF), channels(HD))
    @test same(c.operator, d.operator)
    @test isequal(SQA.simplify(c.rate - d.rate), 0)
  end
end

@testset "the virtual order keeps the Kossakowski matrix positive semidefinite" begin
  C = kossakowski(BF)
  for value in (3.0, 10.0, 100.0)
    numeric = [ComplexF64(SQA.to_complex(SQA.substitute(c, Dict(ω => value)))) for c in C]
    @test eigmin(Hermitian(numeric)) >= -1.0e-12
  end
end

@testset "the virtual order reproduces the exact Floquet exponents through order six" begin
  d, harmonics = 6, 3
  frequencies = [10.0, 14.0, 20.0, 28.0]
  exact = Dict(value => floquet_exponents(d, value, harmonics) for value in frequencies)
  function errors(expansion)
    return [
      slow_distance(exact[value], eigvals(effective_superoperator(expansion, d, value))) for
      value in frequencies
    ]
  end
  before, after = errors(kerr_expansion(6)), errors(BF)
  @test issorted(after; rev=true)
  @test loglog_slope(frequencies, before) >= 5.9
  @test loglog_slope(frequencies, after) >= 6.9
  @test all(after .<= before ./ 5)
end

@testset "lattice invariants vanish on the range of the averaged generator" begin
  data = factorization(BF)
  representation = data.representation
  L0 = data.recurrence.E[1]
  ladder = FE.boson_ladder(representation, L0)
  T = eltype(representation.metric)
  none = zeros(T, size(representation.metric))
  hamiltonian_shift = FE.native_gksl(
    representation, FE.algebra_operator(representation.algebra, [([2, 2], one(T))]), none
  )
  charged_shift = FE.native_gksl(
    representation, FE.algebra_operator(representation.algebra, [([2, 1], one(T))]), none
  )
  for Z in (hamiltonian_shift, charged_shift)
    @test FE.check_invariants(ladder, FE.native_commutator(L0, Z), 8)
  end
  @test_throws FE.VirtualGaugeUnsupported FE.check_invariants(ladder, hamiltonian_shift, 8)
end

@testset "the virtual gauge is refused outside what it certifies" begin
  @test_throws "last retained order 7" kerr_expansion(8)
  for c1 in (-3 // 10, -6 // 5)
    @test_throws "resonant on the Fock lattice" kerr_expansion(7; H=kerr(c1))
  end
  squeezed = kerr(1 // 2, (1 // 20) * (A * A + A' * A'))
  @test_throws "number-diagonal averaged generator" kerr_expansion(7; H=squeezed)
  @test_throws "no static gauge of polynomial degree" kerr_expansion(7; H=squeezed)

  site = FE.boson_site(Rational{Int128})
  algebra = FE.OperatorAlgebra{Complex{Rational{Int128}}}([site, site])
  representation = FE.AlgebraicLiouvilleRepresentation(algebra, 1, 1)
  L0 = FE.AlgebraSuperoperator(
    algebra, Dict{Tuple{Vector{Int},Vector{Int}},eltype(representation.metric)}()
  )
  @test_throws "single bosonic site" FE.boson_ladder(representation, L0)
end

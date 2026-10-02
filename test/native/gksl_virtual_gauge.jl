using Test
using LinearAlgebra: Hermitian, eigmin
using Symbolics: @variables
using FloquetExpansions
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions

@variables ω::Real t::Real

function same(A::Liouvillian, B::Liouvillian)
  return FE.liouvillian_iszero(FE.canonical_liouvillian(FE.SQA.simplify(A - B)))
end
same(A, B) = iszero(FE.SQA.simplify(FE.qadd(A) - FE.qadd(B)))

function kerr_model()
  h = FockSpace(:c)
  a = Destroy(h, :a)
  H = (1 // 2) * a' * a + (3 // 10) * a' * a * a' * a + (2 // 5) * (a + a') * cos(ω * t)
  return H, a
end

function order_six(algorithm)
  H, a = kerr_model()
  return floquet_expansion(
    H, ω, t, GKSLNormalForm(; algorithm), 7; channels=(jump(a, 4 // 5),)
  )
end

# Driven Kerr with loss has no polynomial static gauge at order six: the charged ladder
# solution of [L0, S] = Y never terminates. The last retained order uses a virtual gauge.
const BF = order_six(BlochFeshbach())
const HD = order_six(HoriDeprit())

@testset "order six births the exact newborn channels through the virtual gauge" begin
  _, a = kerr_model()
  for expansion in (BF, HD)
    jumps = channels(expansion)
    @test length(jumps) == 4
    data = getfield(expansion, :completion).factorization
    @test FE.frame_operator(
      data.representation, data.recurrence.channels[2].coefficients[1]
    ) == FE.lower_qadd(data.lowering, FE.qadd(a * a - 2 * a' * a))
    @test isequal(FE.SQA.simplify(jumps[2].rate * ω^4), 72 // 3125)
    @test same(jumps[3].operator, a')
    @test isequal(FE.SQA.simplify(jumps[3].rate * ω^6), 18 // 78125)
    @test same(jumps[4].operator, a' * a' * a * a)
    @test isequal(FE.SQA.simplify(jumps[4].rate * ω^6), 5832 // 78125)
    @test same(
      effective_generator(expansion),
      liouvillian(hamiltonian(expansion); channels=Tuple(jumps)),
    )
    @test getfield(expansion, :completion).virtual_orders == [6]
  end
end

@testset "Bloch-Feshbach and Hori-Deprit agree at the virtual order" begin
  for n in 0:6
    @test same(effective_component(BF, n), effective_component(HD, n))
  end
  @test same(hamiltonian(BF), hamiltonian(HD))
end

@testset "the virtual order keeps the oscillatory micromotion and omits the static part" begin
  for expansion in (BF, HD)
    @test !haskey(micromotion(expansion, 6).components, 0)
    @test !isempty(micromotion(expansion, 6).components)
    @test getfield(expansion, :completion).virtual_orders == [6]
  end
end

@testset "the finite Kossakowski matrix stays positive semidefinite" begin
  C = kossakowski(BF)
  for value in (3.0, 10.0, 100.0)
    numeric = [
      ComplexF64(FE.SQA.to_complex(FE.SQA.substitute(c, Dict(ω => value)))) for c in C
    ]
    @test eigmin(Hermitian(numeric)) >= -1e-12
  end
end

@testset "the lattice invariants of the virtual defect vanish" begin
  data = getfield(BF, :completion).factorization
  rep = data.representation
  ladder = FE.native_boson_ladder(rep, data.recurrence.E[1])
  @test ladder isa FE.BosonLadder
  defect = data.recurrence.virtual[6]
  # E6 - V̂6 is invisible to every lattice invariant ℓ_{m,m'}, well beyond the fitted window.
  @test iszero(FE.max_lattice_invariant(ladder, defect, 24))
  # Positive control: the invariants annihilate every commutator with L0.
  L0 = data.recurrence.E[1]
  H = FE.algebra_operator(rep.algebra, [([2, 2], one(eltype(rep.metric)))])
  Z = FE.native_gksl(rep, H, zeros(eltype(rep.metric), size(rep.metric)))
  @test iszero(FE.max_lattice_invariant(ladder, FE.native_commutator(L0, Z), 8))
  # Negative control: a neutral Hamiltonian shift a†²a² is not invariant-free.
  @test !iszero(FE.max_lattice_invariant(ladder, defect + Z, 8))
  @test !iszero(FE.max_lattice_invariant(ladder, Z, 8))
end

@testset "support detection and the order restriction" begin
  data = getfield(BF, :completion).factorization
  rep = data.representation
  L0 = data.recurrence.E[1]
  @test FE.native_boson_ladder(rep, L0) isa FE.BosonLadder
  T = eltype(rep.metric)
  off_diagonal = FE.native_gksl(
    rep, FE.algebra_operator(rep.algebra, [([2, 1], one(T))]), zeros(T, size(rep.metric))
  )
  @test FE.native_boson_ladder(rep, L0 + off_diagonal) === nothing
  # One order beyond the virtual order: order six needs a static gauge that is not
  # polynomial, but it is no longer the last retained order.
  H, a = kerr_model()
  hd = GKSLNormalForm(; algorithm=HoriDeprit())
  error = try
    floquet_expansion(H, ω, t, hd, 8; channels=(jump(a, 4 // 5),))
    nothing
  catch e
    e
  end
  @test error isa ArgumentError
  @test occursin("explicit static gauge", error.msg)
end

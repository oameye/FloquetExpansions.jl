using Test
using LinearAlgebra: Hermitian, eigmin
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: @variables
using FloquetExpansions
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions
const SQA = SecondQuantizedAlgebra

@variables ω::Real t::Real

function numeric_matrix(matrix, value)
  return [ComplexF64(SQA.to_complex(SQA.substitute(c, Dict(ω => value)))) for c in matrix]
end

function operator_residual(difference::SQA.QAdd, value::Real; include_identity::Bool)
  residual = 0.0
  for (term, coeff) in SQA.simplify(difference)
    include_identity || !isempty(term.ops) || continue
    z = ComplexF64(SQA.to_complex(SQA.substitute(coeff, Dict(ω => value))))
    residual = max(residual, abs(z))
  end
  return residual
end

# The static-gauge-free GKSL normal form is the S = 0 case of graded Gram completion.
# Both realizations are compared in one dissipative frame, where the Kossakowski matrix is the
# complete dissipative datum. The Hamiltonians agree up to a multiple of the identity, which
# the Liouvillian does not see and which the two constructions fix differently.
function agrees_with_gram(H, order, model_channels, algorithm)
  native = floquet_expansion(
    H, ω, t, GKSLNormalForm(; algorithm), order; channels=model_channels
  )
  reference = floquet_expansion(
    H, ω, t, VanVleck(; algorithm), order; channels=model_channels
  )
  frame = dissipative_frame(native)
  completed = positive_completion(reference, Gram(), frame)
  @test dissipative_frame(completed) == frame
  difference = hamiltonian(native) - hamiltonian(completed)
  @test operator_residual(difference, 7.0; include_identity=false) <= 1.0e-10
  @test size(kossakowski(native)) == size(kossakowski(completed))
  for value in (7.0, 25.0, 100.0)
    Cn = numeric_matrix(kossakowski(native), value)
    Cg = numeric_matrix(kossakowski(completed), value)
    @test maximum(abs, Cn - Cg) <= 1.0e-12 * max(1, maximum(abs, Cn))
    @test eigmin(Hermitian(Cn)) >= -1.0e-12
  end
  return native, completed
end

@testset "driven Kerr with loss: GKSLNormalForm equals Gram completion of Van Vleck" begin
  h = FockSpace(:gksl_gram_kerr)
  a = Destroy(h, :a)
  H = (1 // 2) * a' * a + (3 // 10) * a' * a * a' * a + (2 // 5) * (a + a') * cos(ω * t)
  cases = ((HoriDeprit(), (4,)), (BlochFeshbach(), (3,)))
  for (algorithm, orders) in cases, order in orders
    native, completed = agrees_with_gram(H, order, (jump(a, 4 // 5),), algorithm)
    difference = hamiltonian(native) - hamiltonian(completed)
    @test operator_residual(difference, 7.0; include_identity=true) <= 1.0e-10
    @test length(channels(native)) == 1
    @test length(channels(completed)) == 1
    @test only(channels(native)).rate == 4 // 5
    @test size(kossakowski(native)) == (1, 1)
  end
end

@testset "driven qubit with decay: GKSLNormalForm equals Gram completion of Van Vleck" begin
  qubit = PauliSpace(:gksl_gram_qubit)
  σx = Pauli(qubit, :σ, 1)
  σy = Pauli(qubit, :σ, 2)
  σz = Pauli(qubit, :σ, 3)
  σminus = (1 // 2) * (σx - im * σy)
  hd = (HoriDeprit(), (4,))
  models = (
    ((1 // 2) * σz + (2 // 5) * cos(ω * t) * σx, (hd, (BlochFeshbach(), (2,)))),
    (
      (1 // 2) * σz + (1 // 5) * σx + (2 // 5) * cos(ω * t) * σz,
      (hd, (BlochFeshbach(), (3,))),
    ),
  )
  for (H, cases) in models, (algorithm, orders) in cases, order in orders
    agrees_with_gram(H, order, (jump(σminus, 3 // 10),), algorithm)
  end
end

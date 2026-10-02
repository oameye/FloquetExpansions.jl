using Test
using LinearAlgebra: Hermitian, eigmin
using Symbolics: @variables, substitute
using FloquetExpansions
using FloquetExpansions: FloquetExpansions

const FE = FloquetExpansions

@variables ω::Real t::Real

function same(A::Liouvillian, B::Liouvillian)
  return FE.liouvillian_iszero(FE.canonical_liouvillian(FE.SQA.simplify(A - B)))
end
same(A, B) = iszero(FE.SQA.simplify(FE.qadd(A) - FE.qadd(B)))

# The Van Vleck reference runs in Float64 once SQA folds a native factor into a raw
# trigonometric coefficient, so comparisons against it use a tolerance at a sample ω.
small(c) = abs(ComplexF64(FE.SQA.to_complex(FE.SQA.substitute(c, Dict(ω => 10))))) <= 1e-12
function nearly_same(A::Liouvillian, B::Liouvillian)
  difference = FE.canonical_liouvillian(FE.SQA.simplify(A - B))
  return all(small(c) for (_, _, c) in FE.terms(difference))
end
function nearly_same(A, B)
  return all(small(c) for c in values(FE.SQA.simplify(FE.qadd(A) - FE.qadd(B)).arguments))
end

function kerr_hamiltonian()
  h = FockSpace(:c)
  a = Destroy(h, :a)
  H = (1 // 2) * a' * a + (3 // 10) * a' * a * a' * a + (2 // 5) * (a + a') * cos(ω * t)
  return H, a
end

@testset "GKSLNormalForm coincides with VanVleck wherever no static gauge is needed" begin
  H, a = kerr_hamiltonian()
  for algorithm in (HoriDeprit(), BlochFeshbach()), channels in ((), (jump(a, 4 // 5),))
    native = floquet_expansion(H, ω, t, GKSLNormalForm(; algorithm), 4; channels)
    reference = floquet_expansion(H, ω, t, VanVleck(; algorithm), 4; channels)
    for n in 0:3
      @test nearly_same(effective_component(native, n), effective_component(reference, n))
    end
    for n in 1:3
      m, r = micromotion(native, n), micromotion(reference, n)
      for k in union(keys(m.components), keys(r.components))
        @test nearly_same(
          get(m.components, k, zero(m.zero_component)),
          get(r.components, k, zero(r.zero_component)),
        )
      end
    end
  end
end

@testset "GKSLNormalForm driven Kerr is GKSL with a drive-induced channel" begin
  H, a = kerr_hamiltonian()
  native = floquet_expansion(H, ω, t, GKSLNormalForm(), 5; channels=(jump(a, 4 // 5),))
  jumps = channels(native)
  @test length(jumps) == 2
  @test same(jumps[2].operator, a * a - 2 * a' * a)
  @test isequal(FE.SQA.simplify(jumps[2].rate * ω^4), 72 // 3125)
  @test same(
    effective_generator(native), liouvillian(hamiltonian(native); channels=Tuple(jumps))
  )
  # Every finite truncation is completely positive: the Kossakowski matrix is PSD at any ω.
  C = kossakowski(native)
  for value in (3.0, 10.0, 100.0)
    numeric = [
      ComplexF64(FE.SQA.to_complex(FE.SQA.substitute(c, Dict(ω => value)))) for c in C
    ]
    @test eigmin(Hermitian(numeric)) >= -1e-12
  end
  @test_throws ArgumentError positive_completion(native, Gram())
  bf = floquet_expansion(
    H, ω, t, GKSLNormalForm(; algorithm=BlochFeshbach()), 5; channels=(jump(a, 4 // 5),)
  )
  @test all(same(effective_component(bf, n), effective_component(native, n)) for n in 0:4)
end

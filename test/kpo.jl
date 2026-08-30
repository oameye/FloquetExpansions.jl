using Test
using FloquetExpansions
import SecondQuantizedAlgebra as SQA
using Symbolics: Symbolics

# End to end on the spec's running example, the Kerr parametric oscillator: 15 harmonics of a
# quartic Fock-space Hamiltonian, checked against eq:kpo-heff0 and eq:kpo-heff1 as printed.
# Unlike the generic oracle in engine.jl this is a physical system with structure, and it is the
# only test that pins the spec's exact prefactors (17/4, 3/20, 9/40, 20/3, 21/2, 33/4).

hc = FockSpace(:cavity)
a = Destroy(hc, :a)
@variables w::Real δ::Real g₃::Real g₄::Real ξ::Real ξr::Real ξi::Real

# eq:kpo-harmonics, normal ordered, c-number parts dropped: those commute out of both eq:Heff1
# and the conjugation by e^{-iK}, which is why the spec drops them.
function kpo_drive(ξ, ξc, aξ2)
  H0 =
    (δ + 3g₄ * (1 + 2aξ2)) * (a' * a) +
    (3 // 2) * g₄ * (a' * a' * a * a) +
    g₃ * ξ * (a' * a') +
    g₃ * ξc * (a * a)
  H1 =
    g₃ * ((a' * a * a) + (1 + 2aξ2) * a) +
    g₄ * (3ξ * (a' * a' * a) + ξc * (a * a * a) + 3ξ * (1 + aξ2) * a')
  H2 =
    2g₃ * ξ * (a' * a) +
    g₄ * ((a' * a * a * a) + (3 // 2) * (1 + 2aξ2) * (a * a) + (3 // 2) * ξ^2 * (a' * a'))
  H3 = g₃ * ((1 // 3) * (a * a * a) + ξ^2 * a') + 3g₄ * ξ * ((a' * a * a) + (1 + aξ2) * a)
  H4 = g₃ * ξ * (a * a) + g₄ * ((1 // 4) * (a * a * a * a) + 3ξ^2 * (a' * a))
  H5 = g₃ * ξ^2 * a + g₄ * ξ * ((a * a * a) + ξ^2 * a')
  H6 = (3 // 2) * g₄ * ξ^2 * (a * a)
  H7 = g₄ * ξ^3 * a

  comps = Dict{Int,SQA.QAdd}(0 => H0)
  for (m, Hm) in enumerate([H1, H2, H3, H4, H5, H6, H7])
    comps[m] = Hm
    comps[-m] = adjoint(Hm)                    # eq:fourierH
  end
  return PeriodicOperator(comps, w), H0
end

# eq:kpo-heff1, times wd
function kpo_heff1(ξ, aξ2)
  hc_(X) = X + adjoint(X)
  return -((20 // 3) * g₃^2 + 9g₃^2 * aξ2 + (9 // 5) * g₄^2 * (10 + 12aξ2 - 5aξ2^2)) *
         (a' * a) -
         ((10 // 3) * g₃^2 + (9 // 40) * g₄^2 * (85 + 48aξ2)) * (a' * a' * a * a) -
         (17 // 4) * g₄^2 * (a' * a' * a' * a * a * a) - hc_(
    (3 // 20) * g₃ * g₄ * ξ * (105 + 94aξ2) * (a' * a') +
    (21 // 2) * g₃ * g₄ * ξ * (a' * a' * a' * a) +
    (33 // 4) * g₄^2 * ξ^2 * (a' * a' * a' * a'),
  )
end

# "both up to constants" (eq:kpo-heff1): a c-number term is one with an empty operator product.
function dropconst(q::SQA.QAdd)
  out = zero(SQA.QAdd)
  for (term, coeff) in q
    isempty(term.ops) && continue
    out = out + coeff * prod(term.ops)
  end
  return out
end

function maxcoeff(q::SQA.QAdd, subs)
  m = 0.0
  for (_, coeff) in q
    z = SQA.to_num(coeff)
    re = Symbolics.value(Symbolics.substitute(real(z), subs))
    im_ = Symbolics.value(Symbolics.substitute(imag(z), subs))
    m = max(m, abs(complex(Float64(re), Float64(im_))))
  end
  return m
end

@testset "KPO, real pump: eq:kpo-heff0 and eq:kpo-heff1 EXACTLY" begin
  # Real ξ keeps every weight chain clear of the float collapse (D5d), so the spec's exact
  # rationals come back exactly. This is the test that shows they are reproduced, not
  # approximated.
  H, H0 = kpo_drive(ξ, ξ, ξ^2)
  vv = floquet_expansion(H, VanVleck(), 2)

  @test iszero(SQA.simplify(effective_hamiltonian(vv, 0) - H0))          # eq:kpo-heff0
  @test iszero(
    SQA.simplify(dropconst(w * effective_hamiltonian(vv, 1)) - dropconst(kpo_heff1(ξ, ξ^2)))
  )                                                                       # eq:kpo-heff1

  # The spec states the linear terms cancel; nothing else here would notice if they did not.
  got = dropconst(w * effective_hamiltonian(vv, 1))
  @test !any(length(term.ops) == 1 for (term, _) in got)
end

@testset "KPO, complex pump: ξ and ξ* are not interchangeable" begin
  # With ξ real the spec's ξ and ξ* coincide, so a bug swapping them passes unnoticed. Complex
  # ξ separates them. The comparison is numeric here because the larger expression tree pushes
  # a weight chain through SQA's ComplexF64 tier (D5d): the residue is ~1.8e-15, not structural.
  ξc = ξr - im * ξi
  ξz = ξr + im * ξi
  aξ2 = ξr^2 + ξi^2
  H, H0 = kpo_drive(ξz, ξc, aξ2)
  vv = floquet_expansion(H, VanVleck(), 2)

  @test iszero(SQA.simplify(effective_hamiltonian(vv, 0) - H0))
  d = SQA.simplify(
    dropconst(w * effective_hamiltonian(vv, 1)) - dropconst(kpo_heff1(ξz, aξ2))
  )
  @test maxcoeff(d, Dict(g₃ => 0.7, g₄ => 0.4, ξr => 0.6, ξi => 0.9, δ => 0.2)) < 1.0e-12
end

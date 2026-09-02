using FloquetExpansions
using LinearAlgebra: I
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra
using Symbolics: Symbolics

function tomatrix(q::SQA.QAdd, d::Int)
  M = zeros(ComplexF64, d, d)
  for (term, coeff) in q
    z = SQA.to_num(coeff)
    c = complex(Float64(Symbolics.value(real(z))), Float64(Symbolics.value(imag(z))))
    T = Matrix{ComplexF64}(I, d, d)
    for o in term.ops
      E = zeros(ComplexF64, d, d)
      E[o.l1, o.l2] = 1
      T = T * E
    end
    M .+= c .* T
  end
  return M
end

function tomatrix(q::SQA.QAdd, d::Int, substitutions::AbstractDict)
  return tomatrix(SQA.substitute(q, substitutions), d)
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

vanishes(q::SQA.QAdd) = iszero(SQA.simplify(q))
vanishes(X::PeriodicGenerator) = all(vanishes(X[l]) for l in keys(X))

function random_drive(rng, HN, d, M, wd::Symbolics.Num)
  randop() = sum(
    complex(randn(rng), randn(rng)) * SQA.Transition(HN, :σ, i, j) for i in 1:d, j in 1:d
  )
  comps = Dict{Int,SQA.QAdd}()
  h0 = randop()
  comps[0] = h0 + adjoint(h0)
  for m in 1:M
    hm = randop()
    comps[m] = hm
    comps[-m] = adjoint(hm)
  end
  return PeriodicGenerator(comps, wd)
end

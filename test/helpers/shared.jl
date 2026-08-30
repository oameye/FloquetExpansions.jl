using FloquetExpansions
using LinearAlgebra: I
import SecondQuantizedAlgebra as SQA
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
vanishes(X::PeriodicOperator) = all(vanishes(X[l]) for l in keys(X))

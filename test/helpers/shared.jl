using FloquetExpansions
using LinearAlgebra: I, diagm
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra
using Symbolics: Symbolics

function numeric(c::SQA.CNum, substitutions)
  value(x) =
    ComplexF64(Symbolics.value(Symbolics.substitute(x, substitutions; fold=Val(true))))
  z = SQA.to_num(c)
  return value(real(z)) + im * value(imag(z))
end

function operator_matrix(op::SQA.Op, d::Int)
  lowering = diagm(1 => sqrt.(1.0:(d - 1)))
  SQA.is_destroy(op) && return complex(lowering)
  SQA.is_create(op) && return complex(lowering')
  E = zeros(ComplexF64, d, d)
  E[op.l1, op.l2] = 1
  return E
end

function tomatrix(q, d::Int, substitutions::AbstractDict)
  M = zeros(ComplexF64, d, d)
  for (term, c) in FloquetExpansions.qadd(q)
    product = Matrix{ComplexF64}(I, d, d)
    for op in term.ops
      product = product * operator_matrix(op, d)
    end
    M .+= numeric(c, substitutions) .* product
  end
  return M
end

tomatrix(q, d::Int) = tomatrix(q, d, Dict())

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

function same(A::SQA.SuperOperator, B::SQA.SuperOperator)
  return iszero(SQA.simplify(SQA.expand_completeness(A - B)))
end
same(A, B) = iszero(SQA.simplify(FloquetExpansions.qadd(A) - FloquetExpansions.qadd(B)))

channel_map(entry) = SQA.SuperOperator(last(entry) * dissipator(SQA.jump(first(entry))))

function lindbladian(H, channels)
  return foldl(channels; init=hamiltonian_action(H)) do L, (J, rate)
    return L + rate * dissipator(J)
  end
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

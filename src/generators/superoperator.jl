@inline function sandwiches(S::SQA.SuperOperator)
  return (
    (SQA.left_operator(term), SQA.right_operator(term), coefficient) for
    (term, coefficient) in S
  )
end

function canonical_compose(A::SQA.SuperOperator, B::SQA.SuperOperator)
  return SQA.expand_completeness(SQA.compose(A, B))
end

component_commutator(A::SQA.QAdd, B::SQA.QAdd) = SQA.commutator(A, B)
function component_commutator(A::SQA.SuperOperator, B::SQA.SuperOperator)
  return SQA.expand_completeness(SQA.commutator(A, B))
end

"""
    harmonics(L::Lindbladian, ω, t) -> PeriodicGenerator{SuperOperator}
    harmonics(S::SuperOperator, ω, t) -> PeriodicGenerator{SuperOperator}

Decompose a symbolic time-dependent map on density operators into the common
periodic-generator representation. A `Lindbladian` is first lowered to its `SuperOperator`
sandwiches. The decomposition is SecondQuantizedAlgebra's exact `harmonics` of the
`SuperOperator`, with the harmonic label negated to match [`PeriodicGenerator`](@ref). It covers
periodic dependence in a Hamiltonian, jump operator, rate, or any combination.

The components are SecondQuantizedAlgebra `SuperOperator`s with every finite-level
completeness relation applied, so a component is zero exactly when it is the zero map.

See also [`PeriodicGenerator`](@ref), [`floquet_expansion`](@ref).
"""
function harmonics(S::SQA.SuperOperator, w::Symbolics.Num, t::Symbolics.Num)
  components = Dict{Int,SQA.SuperOperator}(
    -n => SQA.expand_completeness(Sn) for (n, Sn) in SQA.harmonics(S, w, t)
  )
  return PeriodicGenerator(components, w)
end

function harmonics(L::SQA.Lindbladian{SQA.QAdd}, w::Symbolics.Num, t::Symbolics.Num)
  return harmonics(SQA.SuperOperator(L), w, t)
end

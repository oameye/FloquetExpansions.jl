@inline function sandwiches(S::SQA.SuperOperator)
  return (
    (SQA.left_operator(term), SQA.right_operator(term), coefficient) for
    (term, coefficient) in S
  )
end

canonical_superoperator(S::SQA.SuperOperator) = SQA.expand_completeness(S)

function canonical_compose(A::SQA.SuperOperator, B::SQA.SuperOperator)
  return canonical_superoperator(SQA.compose(A, B))
end

component_commutator(A::SQA.QAdd, B::SQA.QAdd) = SQA.commutator(A, B)
function component_commutator(A::SQA.SuperOperator, B::SQA.SuperOperator)
  return canonical_superoperator(SQA.commutator(A, B))
end

"""
    harmonics(L::Lindbladian, ω, t) -> PeriodicGenerator{SuperOperator}
    harmonics(S::SuperOperator, ω, t) -> PeriodicGenerator{SuperOperator}

Decompose a symbolic time-dependent map on density operators into the common
periodic-generator representation. A `Lindbladian` is first lowered to its `SuperOperator`
sandwiches. The Fourier decomposition is applied independently to the left and right operator
of every sandwich and to its scalar coefficient, so periodic dependence in a Hamiltonian, jump
operator, rate, or any combination is supported.

The components are SecondQuantizedAlgebra `SuperOperator`s with every finite-level
completeness relation applied, so a component is zero exactly when it is the zero map.

See also [`PeriodicGenerator`](@ref), [`floquet_expansion`](@ref).
"""
function harmonics(S::SQA.SuperOperator, w::Symbolics.Num, t::Symbolics.Num)
  builders = Dict{Int,SQA.SuperBuilder}()
  for (left, right, coefficient) in sandwiches(S)
    left_harmonics = harmonics(left, w, t)
    right_harmonics = harmonics(right, w, t)

    for (left_harmonic, left_component) in component_pairs(left_harmonics),
      (right_harmonic, right_component) in component_pairs(right_harmonics)

      foreach_harmonic_phase(coefficient, w, t) do coefficient_harmonic, phase_coefficient
        harmonic = left_harmonic + right_harmonic + coefficient_harmonic
        builder = get!(SQA.SuperBuilder, builders, harmonic)
        SQA.accumulate!(builder, phase_coefficient, left_component, right_component)
        return nothing
      end
    end
  end
  components = Dict{Int,SQA.SuperOperator}(
    harmonic => canonical_superoperator(SQA.build(builder)) for
    (harmonic, builder) in builders
  )
  return PeriodicGenerator(components, w, zero(SQA.SuperOperator))
end

function harmonics(L::SQA.Lindbladian{SQA.QAdd}, w::Symbolics.Num, t::Symbolics.Num)
  return harmonics(SQA.SuperOperator(L), w, t)
end

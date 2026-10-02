struct NativeExpansionData{T,R,X,O}
  lowering::SQALowering{T,R}
  representation::AlgebraicLiouvilleRepresentation{T,R}
  recurrence::NativeRecurrence{X,O,T}
end

function native_leading_channels(
  representation::AlgebraicLiouvilleRepresentation{T,R}, L0
) where {T,R}
  C0 = native_kossakowski(representation, L0)
  columns, weights = exact_psd_factor(C0, R)
  return Matrix{T}(columns), T.(weights)
end

function native_input_degree(lowering::SQALowering, L::Dict)
  degree = 0
  for (_, X) in L, ((left, right), _) in X.terms
    degree = max(
      degree,
      monomial_degree(lowering.algebra, left),
      monomial_degree(lowering.algebra, right),
    )
  end
  return degree
end

function native_expansion_data(
  ::Type{T},
  algorithm::Union{BlochFeshbach,HoriDeprit},
  generator::PeriodicGenerator{Liouvillian},
  N::Int,
  parameters::AbstractDict,
  gauge_degree::Int,
) where {T}
  lowering = SQALowering{T}(generator, parameters)
  L = lower_generator(lowering, generator)
  haskey(L, 0) || (L[0] = zero(first(values(L))))
  degree = native_input_degree(lowering, L) + N
  while true
    representation = AlgebraicLiouvilleRepresentation(
      lowering.algebra, degree, gauge_degree
    )
    leading, weights = native_leading_channels(representation, L[0])
    try
      recurrence = native_recurrence(
        algorithm, representation, NoHomologicalInverse(), L, leading, weights, N, 0.0
      )
      return NativeExpansionData(lowering, representation, recurrence)
    catch error
      error isa NativeFrameError || rethrow()
      degree += 2
    end
  end
end

function native_expansion_data(
  algorithm::Union{BlochFeshbach,HoriDeprit},
  generator::PeriodicGenerator{Liouvillian},
  N::Int,
  parameters::AbstractDict,
  gauge_degree::Int,
)
  try
    return native_expansion_data(
      Complex{Rational{Int128}}, algorithm, generator, N, parameters, gauge_degree
    )
  catch error
    error isa OverflowError || rethrow()
    return native_expansion_data(
      Complex{Rational{BigInt}}, algorithm, generator, N, parameters, gauge_degree
    )
  end
end

function lift_superoperator(lowering::SQALowering, S::AlgebraSuperoperator)
  result = zero(Liouvillian)
  for ((X, Y), c) in S.terms
    left = lift_operator(lowering, algebra_operator(lowering.algebra, [(X, 1)]))
    right = lift_operator(lowering, algebra_operator(lowering.algebra, [(Y, 1)]))
    result = result + sqa_scalar(c) * action(left, right)
  end
  return result
end

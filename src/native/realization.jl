const NATIVE_GAUGE_DEGREE = 3

struct NativeRealization{F<:DissipativeFrame,J<:RateWeightedJump,X} <: Completion
  frame::F
  retained_kossakowski::Vector{KossakowskiMatrix}
  kossakowski::KossakowskiMatrix
  channels::Vector{J}
  positivity_conditions::Vector{SQA.CNum}
  regularity_conditions::Vector{SQA.CNum}
  factorization::X
  hamiltonian::SQA.QAdd
  generator::Liouvillian
end

function lift_harmonics(lift::CoefficientLift, harmonics::Dict, wd::Symbolics.Num)
  components = Dict{Int,Liouvillian}(
    k => lift_superoperator(lift, X) for (k, X) in harmonics if !iszero(X)
  )
  return PeriodicGenerator(components, wd, zero(Liouvillian))
end

function native_micromotion_harmonics(data::NativeExpansionData, ::HoriDeprit, N::Int)
  return [data.recurrence.kick[n + 1] for n in 1:N]
end

function native_micromotion_harmonics(data::NativeExpansionData, ::BlochFeshbach, N::Int)
  L0 = data.recurrence.E[1]
  generator = kick_log(data.recurrence.kick, one(L0), N)
  return [generator[n + 1] for n in 1:N]
end

function native_micromotion(
  data::NativeExpansionData, lift::CoefficientLift, algorithm, wd, N::Int
)
  return [
    lift_harmonics(lift, harmonics, wd) for
    harmonics in native_micromotion_harmonics(data, algorithm, N)
  ]
end

function native_channel_operator(
  data::NativeExpansionData, lift::CoefficientLift, channel::GradedChannel, wd, N::Int
)
  operator = zero(SQA.QAdd)
  for (k, coefficients) in pairs(channel.coefficients)
    order = k - 1
    order + channel.onset > N && break
    lifted = lift_operator(lift, frame_operator(data.representation, coefficients))
    operator = operator + reattach(lifted, wd, order)
  end
  return SQA.simplify(operator)
end

function native_channels(data::NativeExpansionData, lift::CoefficientLift, wd, N::Int)
  return RateWeightedJump{SQA.QAdd}[
    jump(
      native_channel_operator(data, lift, channel, wd, N),
      lifted_coefficient(lift, channel.weight) *
      convert(SQA.CNum, inverse_drive_power(wd, channel.onset)),
    ) for channel in data.recurrence.channels if channel.onset <= N
  ]
end

function native_support(data::NativeExpansionData)
  support = Int[]
  for channel in data.recurrence.channels, coefficients in channel.coefficients
    for (k, c) in pairs(coefficients)
      iszero(c) || k in support || push!(support, k)
    end
  end
  return sort!(support)
end

function native_frame(data::NativeExpansionData, support::Vector{Int})
  operators = [
    lift_operator(
      data.lowering,
      algebra_operator(data.lowering.algebra, [(data.representation.frame[k], 1)]),
    ) for k in support
  ]
  return DissipativeFrame(operators)
end

function retained_amplitudes(channel::GradedChannel, N::Int)
  return [
    (k - 1, c) for (k, c) in pairs(channel.coefficients) if k - 1 + channel.onset <= N
  ]
end

function native_finite_kossakowski(data::NativeExpansionData{T}, N::Int) where {T}
  n = length(data.representation.frame)
  graded = Dict{Int,Matrix{T}}()
  for channel in data.recurrence.channels
    channel.onset <= N || continue
    amplitudes = retained_amplitudes(channel, N)
    for (k, l) in amplitudes, (q, r) in amplitudes
      C = get!(() -> exact_zeros(T, n, n), graded, channel.onset + k + q)
      add_outer!(C, channel.weight, l, r)
    end
  end
  return graded
end

function graded_coefficient(entries, wd::Symbolics.Num, coefficient)
  result = convert(SQA.CNum, 0)
  for (order, c) in entries
    iszero(c) && continue
    scale = convert(SQA.CNum, inverse_drive_power(wd, order))
    result = result + convert(SQA.CNum, coefficient(c)) * scale
  end
  return result
end

function support_matrix(graded, support::Vector{Int}, wd::Symbolics.Num, coefficient)
  m = length(support)
  K = KossakowskiMatrix(undef, m, m)
  for j in 1:m, i in 1:m
    K[i, j] = graded_coefficient(
      ((order, C[support[i], support[j]]) for (order, C) in graded), wd, coefficient
    )
  end
  return K
end

function check_support(C::AbstractMatrix, inside::AbstractVector{Bool})
  for index in CartesianIndices(C)
    μ, ν = Tuple(index)
    iszero(C[index]) ||
      (inside[μ] && inside[ν]) ||
      throw(ArgumentError("a retained Kossakowski component leaves the channel frame"))
  end
  return C
end

function retained_native_kossakowski(
  data::NativeExpansionData, support::Vector{Int}, wd::Symbolics.Num, N::Int, coefficient
)
  inside = falses(length(data.representation.frame))
  inside[support] .= true
  return KossakowskiMatrix[
    support_matrix(
      Dict(
        n => check_support(
          native_kossakowski(data.representation, data.recurrence.E[n + 1]), inside
        ),
      ),
      support,
      wd,
      coefficient,
    ) for n in 0:N
  ]
end

function native_finite_superoperators(
  data::NativeExpansionData{T,R}, graded, N::Int
) where {T,R}
  representation = data.representation
  zero_frame = exact_zeros(T, length(representation.frame), length(representation.frame))
  zero_hamiltonian = AlgebraOperator{T,R}(representation.algebra, Dict{Monomial,T}())
  return Dict(
    order => native_gksl(
      representation,
      order <= N ? data.recurrence.H[order + 1] : zero_hamiltonian,
      get(graded, order, zero_frame),
    ) for order in union(0:N, keys(graded))
  )
end

function native_finite_generator(
  data::NativeExpansionData, lift::CoefficientLift, graded, wd, N::Int
)
  generator = Liouvillian(LiouvillianTerms())
  for (order, G) in native_finite_superoperators(data, graded, N)
    scale = convert(SQA.CNum, inverse_drive_power(wd, order))
    lift_superoperator!(generator, lift, G, scale)
  end
  return generator
end

function floquet_expansion_impl(
  generator::PeriodicGenerator{Liouvillian},
  gauge::GKSLNormalForm,
  order::Int,
  provenance::R,
) where {R<:FloquetProvenance}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  data = try
    native_expansion_data(gauge.algorithm, generator, order - 1, NATIVE_GAUGE_DEGREE)
  catch error
    error isa NativePositivityError || rethrow()
    throw(
      ArgumentError(
        "order $order of the GKSL normal form has no static gauge of polynomial degree " *
        "at most $NATIVE_GAUGE_DEGREE; lower the requested order",
      ),
    )
  end
  return native_floquet_expansion(data, generator, gauge, provenance)
end

native_positivity_conditions(::ExactCoefficients) = SQA.CNum[]

native_factorization(data::NativeExpansionData, ::ExactCoefficients) = data

function native_hamiltonian_series(data::NativeExpansionData, lift::CoefficientLift, wd)
  H = zero(SQA.QAdd)
  for (n, Hn) in pairs(data.recurrence.H)
    H = H + reattach(lift_operator(lift, Hn), wd, n - 1)
  end
  return SQA.simplify(H)
end

function native_floquet_expansion(
  data::NativeExpansionData,
  generator::PeriodicGenerator{Liouvillian},
  gauge::GKSLNormalForm,
  provenance::R,
) where {R<:FloquetProvenance}
  return native_floquet_expansion(data, generator, gauge, provenance, ExactCoefficients())
end

function native_floquet_expansion(
  data::NativeExpansionData,
  generator::PeriodicGenerator{Liouvillian},
  gauge::GKSLNormalForm,
  provenance::R,
  coefficients,
) where {R<:FloquetProvenance}
  order = length(data.recurrence.E)
  N = order - 1
  wd = generator.wd
  lift = CoefficientLift(data.lowering, coefficients)
  effective = Liouvillian[lift_superoperator(lift, E) for E in data.recurrence.E]
  micromotion = native_micromotion(data, lift, gauge.algorithm, wd, N)
  H = native_hamiltonian_series(data, lift, wd)
  jumps = native_channels(data, lift, wd, N)
  support = native_support(data)
  frame = native_frame(data, support)
  graded = native_finite_kossakowski(data, N)
  finite = native_finite_generator(data, lift, graded, wd, N)
  retained = retained_native_kossakowski(data, support, wd, N, coefficients)
  realization = NativeRealization(
    frame,
    retained,
    support_matrix(graded, support, wd, coefficients),
    jumps,
    native_positivity_conditions(coefficients),
    SQA.CNum[],
    native_factorization(data, coefficients),
    H,
    finite,
  )
  return FloquetExpansion(
    generator, micromotion, effective, gauge, order, realization, provenance
  )
end

function floquet_expansion_impl(
  generator::PeriodicGenerator{SQA.QAdd}, gauge::GKSLNormalForm, order::Int, provenance::R
) where {R<:FloquetProvenance}
  reference = floquet_expansion_impl(
    generator, VanVleck(gauge.algorithm), order, provenance
  )
  return FloquetExpansion(
    generator,
    getfield(reference, :micromotion_components),
    getfield(reference, :effective_components),
    gauge,
    order,
    Uncompleted(),
    provenance,
  )
end

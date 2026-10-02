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

function lift_harmonics(lowering::SQALowering, harmonics::Dict, wd::Symbolics.Num)
  components = Dict{Int,Liouvillian}(
    k => lift_superoperator(lowering, X) for (k, X) in harmonics if !iszero(X)
  )
  return PeriodicGenerator(components, wd, zero(Liouvillian))
end

function native_micromotion(data::NativeExpansionData, ::HoriDeprit, wd, N::Int)
  return [lift_harmonics(data.lowering, data.recurrence.kick[n + 1], wd) for n in 1:N]
end

function native_micromotion(data::NativeExpansionData, ::BlochFeshbach, wd, N::Int)
  L0 = data.recurrence.E[1]
  generator = kick_log(data.recurrence.kick, one(L0), N)
  return [lift_harmonics(data.lowering, generator[n + 1], wd) for n in 1:N]
end

function native_channel_operator(
  data::NativeExpansionData, channel::GradedChannel, wd, N::Int
)
  operator = zero(SQA.QAdd)
  for (k, coefficients) in pairs(channel.coefficients)
    order = k - 1
    order + channel.onset > N && break
    lifted = lift_operator(data.lowering, frame_operator(data.representation, coefficients))
    operator = operator + reattach(lifted, wd, order)
  end
  return SQA.simplify(operator)
end

function native_channels(data::NativeExpansionData, wd, N::Int)
  return RateWeightedJump{SQA.QAdd}[
    jump(
      native_channel_operator(data, channel, wd, N),
      sqa_scalar(channel.weight) * inverse_drive_power(wd, channel.onset),
    ) for channel in data.recurrence.channels if channel.onset <= N
  ]
end

function native_frame(data::NativeExpansionData)
  support = Int[]
  for channel in data.recurrence.channels, coefficients in channel.coefficients
    for (k, c) in pairs(coefficients)
      iszero(c) || k in support || push!(support, k)
    end
  end
  sort!(support)
  operators = [
    lift_operator(
      data.lowering,
      algebra_operator(data.lowering.algebra, [(data.representation.frame[k], 1)]),
    ) for k in support
  ]
  return DissipativeFrame(operators)
end

function floquet_expansion_impl(
  generator::PeriodicGenerator{Liouvillian},
  gauge::GKSLNormalForm,
  order::Int,
  provenance::R,
) where {R<:FloquetProvenance}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  N = order - 1
  wd = generator.wd
  data = native_expansion_data(gauge.algorithm, generator, N, Dict{Any,Any}(), 3)
  effective = Liouvillian[lift_superoperator(data.lowering, E) for E in data.recurrence.E]
  micromotion = native_micromotion(data, gauge.algorithm, wd, N)
  H = zero(SQA.QAdd)
  for (n, Hn) in pairs(data.recurrence.H)
    H = H + reattach(lift_operator(data.lowering, Hn), wd, n - 1)
  end
  H = SQA.simplify(H)
  jumps = native_channels(data, wd, N)
  finite = liouvillian(H; channels=Tuple(jumps))
  frame = native_frame(data)
  retained = KossakowskiMatrix[
    kossakowski(reattach(effective[n + 1], wd, n), frame) for n in 0:N
  ]
  realization = NativeRealization(
    frame,
    retained,
    kossakowski(finite, frame),
    jumps,
    SQA.CNum[],
    SQA.CNum[],
    data,
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

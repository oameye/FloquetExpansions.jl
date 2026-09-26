struct JumpAmplitudeSeed
  amplitude::PeriodicGenerator{SQA.QAdd}
  rate::SQA.CNum
  reference::DissipativeSeedRef
end

struct TransportedJumpAmplitude
  seed::JumpAmplitudeSeed
  coefficients::Vector{PeriodicGenerator{SQA.QAdd}}
end

struct HarmonicJumpChannel
  operator::SQA.QAdd
  rate::SQA.CNum
  reference::DissipativeSeedRef
  harmonic::Int
end

struct NativeGKSLExpansion{F<:FloquetExpansion}
  coherent::F
  amplitudes::Vector{TransportedJumpAmplitude}
  channels::Vector{HarmonicJumpChannel}
  generator::Liouvillian
end

function static_jump_rate(
  rate::SQA.CNum, operator::SQA.QAdd, wd::Symbolics.Num, t::Symbolics.Num
)::SQA.CNum
  rate_harmonics = harmonics((rate * one(operator))::SQA.QAdd, wd, t)
  any(!iszero(harmonic) for harmonic in keys(rate_harmonics)) && throw(
    ArgumentError(
      "the native GKSL expansion requires a time-independent jump rate; " *
      "pass a periodic `collapse` operator that carries the square root of the rate instead",
    ),
  )
  return rate
end

function jump_amplitude_seed(
  provenance::MicroscopicProvenance,
  reference::DissipativeSeedRef,
  wd::Symbolics.Num,
  t::Symbolics.Num,
)
  if reference.kind == COLLAPSE_SEED
    operator = provenance.collapse_operators[reference.index]
    return JumpAmplitudeSeed(harmonics(operator, wd, t), convert(SQA.CNum, 1), reference)
  end
  operator = provenance.jump_operators[reference.index]
  rate = static_jump_rate(provenance.jump_rates[reference.index], operator, wd, t)
  return JumpAmplitudeSeed(harmonics(operator, wd, t), rate, reference)
end

function jump_amplitude_seeds(
  channels::LiouvillianChannelCollection, wd::Symbolics.Num, t::Symbolics.Num
)
  provenance = microscopic_provenance(channels)
  result = JumpAmplitudeSeed[]
  sizehint!(result, length(provenance.order))
  for reference in provenance.order
    push!(result, jump_amplitude_seed(provenance, reference, wd, t))
  end
  return result
end

function transport_jump_amplitude(
  seed::JumpAmplitudeSeed,
  micromotion_components::Vector{PeriodicGenerator{SQA.QAdd}},
  order::Int,
)
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  length(micromotion_components) == order - 1 || throw(
    ArgumentError("the coherent micromotion truncation is inconsistent with the order")
  )

  amplitude = seed.amplitude
  nodes = (order * (order + 1)) ÷ 2
  dressed = [zero(amplitude) for _ in 1:nodes]
  coefficients = PeriodicGenerator{SQA.QAdd}[]
  sizehint!(coefficients, order)

  for n in 0:(order - 1)
    dressed[triindex(n, 0)] = iszero(n) ? amplitude : zero(amplitude)
    for j in 1:n
      dressed[triindex(n, j)] = SQA.simplify(
        dressed_generator_node(micromotion_components, dressed, n, j, amplitude)
      )::PeriodicGenerator{SQA.QAdd}
    end

    coefficient = zero(amplitude)
    for j in 0:n
      coefficient = coefficient + weight_generator(amplitude, j) * dressed[triindex(n, j)]
    end
    push!(coefficients, SQA.simplify(coefficient)::PeriodicGenerator{SQA.QAdd})
  end

  return TransportedJumpAmplitude(seed, coefficients)
end

function finite_transported_amplitude(amplitude::TransportedJumpAmplitude)
  result = zero(amplitude.seed.amplitude)
  for (index, coefficient) in enumerate(amplitude.coefficients)
    result = result + reattach(coefficient, index - 1)
  end
  return SQA.simplify(result)::PeriodicGenerator{SQA.QAdd}
end

function harmonic_jump_channels(amplitudes::Vector{TransportedJumpAmplitude})
  result = HarmonicJumpChannel[]
  for amplitude in amplitudes
    seed = amplitude.seed
    iszero(seed.rate) && continue
    finite = finite_transported_amplitude(amplitude)
    harmonic_indices = sort!(collect(keys(finite)))
    sizehint!(result, length(result) + length(harmonic_indices))
    for harmonic in harmonic_indices
      operator = SQA.simplify(finite[harmonic])::SQA.QAdd
      iszero(operator) && continue
      push!(result, HarmonicJumpChannel(operator, seed.rate, seed.reference, harmonic))
    end
  end
  return result
end

@inline channel_liouvillian(channel::HarmonicJumpChannel) =
  channel.rate * dissipator(channel.operator)

function native_gksl_generator(
  coherent::FloquetExpansion{G,P,SQA.QAdd}, channels::Vector{HarmonicJumpChannel}
) where {G,P<:PeriodicGenerator{SQA.QAdd}}
  generator = hamiltonian_action(effective_generator(coherent)::SQA.QAdd)
  for channel in channels
    generator = generator + channel_liouvillian(channel)
  end
  return SQA.simplify(generator)::Liouvillian
end

function native_dissipative_component(
  amplitudes::Vector{TransportedJumpAmplitude}, grade::Int
)::Liouvillian
  grade >= 0 || throw(ArgumentError("grade must be nonnegative"))
  result = zero(Liouvillian)

  for amplitude in amplitudes
    retained = length(amplitude.coefficients)
    grade < 2 * retained - 1 || continue
    rate = amplitude.seed.rate
    iszero(rate) && continue
    for left_grade in 0:grade
      right_grade = grade - left_grade
      left_grade < retained || continue
      right_grade < retained || continue
      left = amplitude.coefficients[left_grade + 1]
      right = amplitude.coefficients[right_grade + 1]
      for harmonic in keys(left)
        right_component = right[harmonic]
        iszero(right_component) && continue
        result = result + rate * cross_dissipator(left[harmonic], right_component)
      end
    end
  end

  return SQA.simplify(result)::Liouvillian
end

function native_gksl_expansion(
  H::SQA.QField,
  wd::Symbolics.Num,
  t::Symbolics.Num,
  order::Int,
  channels::LiouvillianChannelCollection,
)
  coherent = floquet_expansion(harmonics(qadd(H), wd, t), VanVleck(), order)
  micromotion_components = getfield(coherent, :micromotion_components)

  seeds = jump_amplitude_seeds(channels, wd, t)
  amplitudes = TransportedJumpAmplitude[]
  sizehint!(amplitudes, length(seeds))
  for seed in seeds
    push!(amplitudes, transport_jump_amplitude(seed, micromotion_components, order))
  end

  harmonic_channels = harmonic_jump_channels(amplitudes)
  generator = native_gksl_generator(coherent, harmonic_channels)
  return NativeGKSLExpansion(coherent, amplitudes, harmonic_channels, generator)
end

const CPPeriodicAmplitude = PeriodicGenerator{SQA.QAdd}

# Physical amplitude data are kept separate from generic Liouvillian algebra. The rate is a
# nonnegative physical channel weight; all drive dependence belongs to `amplitude` in this first
# implementation tranche.
struct PhysicalAmplitudeSeed
  amplitude::CPPeriodicAmplitude
  rate::SQA.CNum
  kind::DissipativeSeedKind
  source_index::Int
end

struct TransportedAmplitudeSeries
  seed::PhysicalAmplitudeSeed
  coefficients::Vector{CPPeriodicAmplitude}
end

# One finite Kraus row after coherent-frame transport and period averaging. Keeping the physical
# seed and Fourier harmonic avoids lowering provenance into an opaque dense Kossakowski matrix.
struct CPAmplitudeChannel
  operator::SQA.QAdd
  rate::SQA.CNum
  kind::DissipativeSeedKind
  source_index::Int
  harmonic::Int
end

struct CPHFEReconstruction{F}
  coherent::F
  amplitudes::Vector{TransportedAmplitudeSeries}
  channels::Vector{CPAmplitudeChannel}
  generator::Liouvillian
end

function static_jump_rate(
  channel::RateWeightedJump, wd::Symbolics.Num, t::Symbolics.Num
)::SQA.CNum
  prototype = one(qadd(channel.operator))
  rate_harmonics = harmonics(channel.rate * prototype, wd, t)
  any(!iszero(harmonic) for harmonic in keys(rate_harmonics)) && throw(
    ArgumentError(
      "native CP-HFE amplitude transport currently requires time-independent jump rates; " *
      "put periodic amplitude dependence in the jump operator instead",
    ),
  )
  return channel.rate
end

function physical_amplitude_seed(
  channel::CollapseChannel, source_index::Int, wd::Symbolics.Num, t::Symbolics.Num
)
  amplitude = harmonics(qadd(channel.operator), wd, t)
  return PhysicalAmplitudeSeed(
    amplitude, convert(SQA.CNum, 1), COLLAPSE_SEED, source_index
  )
end

function physical_amplitude_seed(
  channel::RateWeightedJump, source_index::Int, wd::Symbolics.Num, t::Symbolics.Num
)
  rate = static_jump_rate(channel, wd, t)
  amplitude = harmonics(qadd(channel.operator), wd, t)
  return PhysicalAmplitudeSeed(amplitude, rate, JUMP_SEED, source_index)
end

function physical_amplitude_seed(
  channel, source_index::Int, wd::Symbolics.Num, t::Symbolics.Num
)
  channel isa LiouvillianChannel || throw(
    ArgumentError("channels must contain only `collapse(...)` and `jump(...)` values")
  )
  return physical_amplitude_seed(channel, source_index, wd, t)
end

function physical_amplitude_seeds(
  channels::LiouvillianChannelCollection, wd::Symbolics.Num, t::Symbolics.Num
)
  result = PhysicalAmplitudeSeed[]
  sizehint!(result, length(channels))
  for (source_index, channel) in enumerate(channels)
    push!(result, physical_amplitude_seed(channel, source_index, wd, t))
  end
  return result
end

function transport_amplitude_series(
  seed::PhysicalAmplitudeSeed, kicks::Vector{CPPeriodicAmplitude}, order::Int
)
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  length(kicks) == order - 1 ||
    throw(ArgumentError("coherent kick truncation is inconsistent with requested order"))

  amplitude = seed.amplitude
  nodes = (order * (order + 1)) ÷ 2
  dressed = [zero(amplitude) for _ in 1:nodes]
  coefficients = CPPeriodicAmplitude[]
  sizehint!(coefficients, order)

  for n in 0:(order - 1)
    dressed[triindex(n, 0)] = iszero(n) ? amplitude : zero(amplitude)
    for j in 1:n
      dressed[triindex(n, j)] = SQA.simplify(
        dressed_generator_node(kicks, dressed, n, j, amplitude)
      )::CPPeriodicAmplitude
    end

    coefficient = zero(amplitude)
    for j in 0:n
      coefficient = coefficient + weight_generator(amplitude, j) * dressed[triindex(n, j)]
    end
    push!(coefficients, SQA.simplify(coefficient)::CPPeriodicAmplitude)
  end

  return TransportedAmplitudeSeries(seed, coefficients)
end

function finite_transported_amplitude(series::TransportedAmplitudeSeries)
  amplitude = zero(series.seed.amplitude)
  for (index, coefficient) in enumerate(series.coefficients)
    amplitude = amplitude + reattach(coefficient, index - 1)
  end
  return SQA.simplify(amplitude)::CPPeriodicAmplitude
end

function reconstruct_cp_amplitude_channels(amplitudes::Vector{TransportedAmplitudeSeries})
  result = CPAmplitudeChannel[]
  for series in amplitudes
    finite = finite_transported_amplitude(series)
    harmonic_indices = sort!(collect(keys(finite)))
    sizehint!(result, length(result) + length(harmonic_indices))
    for harmonic in harmonic_indices
      operator = SQA.simplify(finite[harmonic])::SQA.QAdd
      iszero(operator) && continue
      seed = series.seed
      push!(
        result,
        CPAmplitudeChannel(
          operator, seed.rate, seed.kind, seed.source_index, harmonic
        ),
      )
    end
  end
  return result
end

@inline function cp_amplitude_channel_liouvillian(channel::CPAmplitudeChannel)
  return channel.rate * dissipator(channel.operator)
end

function reconstruct_cp_effective_generator(
  coherent::FloquetExpansion, channels::Vector{CPAmplitudeChannel}
)
  H_eff = effective_generator(coherent)::SQA.QAdd
  generator = hamiltonian_action(H_eff)
  for channel in channels
    generator = generator + cp_amplitude_channel_liouvillian(channel)
  end
  return SQA.simplify(generator)::Liouvillian
end

# Retained coefficient of the period-averaged dissipative Gram/Kraus form before the finite square
# is evaluated. This is diagnostic/reference data; the finite CP generator is built from complete
# retained amplitudes and is deliberately not truncated after squaring.
function cp_dissipative_component(
  amplitudes::Vector{TransportedAmplitudeSeries}, grade::Int
)::Liouvillian
  grade >= 0 || throw(ArgumentError("grade must be nonnegative"))
  result = zero(Liouvillian)

  for series in amplitudes
    grade < 2 * length(series.coefficients) - 1 || continue
    rate = series.seed.rate
    for left_grade in 0:grade
      right_grade = grade - left_grade
      left_grade < length(series.coefficients) || continue
      right_grade < length(series.coefficients) || continue
      left = series.coefficients[left_grade + 1]
      right = series.coefficients[right_grade + 1]
      for harmonic in keys(left)
        right_component = right[harmonic]
        iszero(right_component) && continue
        result =
          result + rate * cross_dissipator(left[harmonic], right_component)
      end
    end
  end

  return SQA.simplify(result)::Liouvillian
end

function cp_hfe_reconstruction(
  H::SQA.QField,
  wd::Symbolics.Num,
  t::Symbolics.Num,
  order::Int,
  channels::LiouvillianChannelCollection,
)
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  coherent_generator = harmonics(qadd(H), wd, t)
  require_hermitian_drive(coherent_generator)
  coherent = floquet_expansion(coherent_generator, VanVleck(), order)
  kicks = getfield(coherent, :kick_components)::Vector{CPPeriodicAmplitude}

  seeds = physical_amplitude_seeds(channels, wd, t)
  amplitudes = TransportedAmplitudeSeries[]
  sizehint!(amplitudes, length(seeds))
  for seed in seeds
    push!(amplitudes, transport_amplitude_series(seed, kicks, order))
  end

  reconstructed_channels = reconstruct_cp_amplitude_channels(amplitudes)
  generator = reconstruct_cp_effective_generator(coherent, reconstructed_channels)
  return CPHFEReconstruction(coherent, amplitudes, reconstructed_channels, generator)
end

if !isdefined(@__MODULE__, :GramStinespringWord)
  include(joinpath(@__DIR__, "gram_stinespring_triangle_reference.jl"))
end

function physical_gram_metric_series(
  words::Vector{GramStinespringWord{T}}, generator_order::Int, zero_component::T
) where {T}
  generator_order >= 0 || throw(ArgumentError("generator order must be nonnegative"))
  max_channel_order = 2 * (generator_order + 1)
  metric = [formal_period_zero(zero_component) for _ in 0:max_channel_order]

  for left in words, right in words
    left.output_number == right.output_number || continue
    channel_order = left.output_number + left.drift_order + right.drift_order
    overlap_right_left = mixed_qr_gram_polynomial(left.metadata, right.metadata)
    all(iszero, overlap_right_left) && continue
    system_term = adjoint(left.system_value) * right.system_value
    formal_period_accumulate!(
      metric[channel_order + 1], conj.(overlap_right_left), system_term
    )
  end
  return metric
end

function physical_gram_channel_series(
  words::Vector{GramStinespringWord{T}}, generator_order::Int, zero_superoperator::S; paste
) where {T,S}
  max_channel_order = 2 * (generator_order + 1)
  channel = [formal_period_zero(zero_superoperator) for _ in 0:max_channel_order]

  for left in words, right in words
    left.output_number == right.output_number || continue
    channel_order = left.output_number + left.drift_order + right.drift_order
    overlap_right_left = mixed_qr_gram_polynomial(left.metadata, right.metadata)
    all(iszero, overlap_right_left) && continue
    system_term = paste(left.system_value, right.system_value)
    formal_period_accumulate!(channel[channel_order + 1], overlap_right_left, system_term)
  end
  return channel
end

function physical_word_harmonic(word::GramStinespringWord)
  return sum(vertex.harmonic for vertex in word.metadata; init=0)
end

function physical_gram_channel_order_by_phase(
  words::Vector{GramStinespringWord{T}},
  requested_order::Int,
  zero_superoperator::S;
  paste,
) where {T,S}
  requested_order >= 0 || throw(ArgumentError("requested order must be nonnegative"))
  resolved = Dict{Int,FormalPeriodMatrix{S}}()

  for left in words, right in words
    left.output_number == right.output_number || continue
    channel_order = left.output_number + left.drift_order + right.drift_order
    channel_order == requested_order || continue
    overlap_right_left = mixed_qr_gram_polynomial(left.metadata, right.metadata)
    all(iszero, overlap_right_left) && continue

    phase_harmonic = physical_word_harmonic(left) - physical_word_harmonic(right)
    component = get!(resolved, phase_harmonic) do
      formal_period_zero(zero_superoperator)
    end
    system_term = paste(left.system_value, right.system_value)
    formal_period_accumulate!(component, overlap_right_left, system_term)
  end
  return resolved
end

function physical_log_second_order_phase_average(
  first_order::FormalPeriodMatrix{T},
  second_order_zero_phase::FormalPeriodMatrix{T},
  zero_component::T,
) where {T}
  square = formal_period_product(first_order, first_order, zero_component)
  result = formal_period_zero(zero_component)
  formal_period_add!(result, second_order_zero_phase, zero_component)
  formal_period_add!(result, formal_period_scale(-1 // 2, square, zero_component), zero_component)
  return result
end

function physical_rr_by_first_mismatch(
  words::Vector{GramStinespringWord{T}}, zero_superoperator::S; paste
) where {T,S}
  resolved = Dict{Int,FormalPeriodMatrix{S}}()

  for left in words, right in words
    left.output_number == 2 || continue
    right.output_number == 2 || continue
    iszero(left.drift_order) || continue
    iszero(right.drift_order) || continue
    physical_word_harmonic(left) == physical_word_harmonic(right) || continue

    overlap_right_left = mixed_qr_gram_polynomial(left.metadata, right.metadata)
    all(iszero, overlap_right_left) && continue
    mismatch = abs(left.metadata[1].harmonic - right.metadata[1].harmonic)
    component = get!(resolved, mismatch) do
      formal_period_zero(zero_superoperator)
    end
    system_term = paste(left.system_value, right.system_value)
    formal_period_accumulate!(component, overlap_right_left, system_term)
  end
  return resolved
end

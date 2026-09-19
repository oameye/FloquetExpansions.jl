@enum SimplexSlowClass begin
  SimplexSlowPrimitive
  SimplexSlowReducible
  SimplexOutsideSlowSector
  SimplexChannelOrthogonal
end

struct SimplexOutputWord{C,H,N}
  channels::NTuple{N,C}
  harmonics::NTuple{N,H}
end

struct SimplexSlowOverlap{T}
  class::SimplexSlowClass
  coefficient::T
end

simplex_output_number(::SimplexOutputWord{C,H,N}) where {C,H,N} = N

function simplex_total_sideband(word::SimplexOutputWord{C,H,N}) where {C,H,N}
  total = zero(H)
  for harmonic in word.harmonics
    total += harmonic
  end
  return total
end

function simplex_slow_overlap(
  left::SimplexOutputWord{C,H,N},
  right::SimplexOutputWord{C,H,N};
  imaginary,
  inverse_weight,
) where {C,H,N}
  zero_coefficient = zero(imaginary)
  left.channels == right.channels ||
    return SimplexSlowOverlap(SimplexChannelOrthogonal, zero_coefficient)

  difference = ntuple(index -> left.harmonics[index] - right.harmonics[index], N)
  total = zero(H)
  for harmonic in difference
    total += harmonic
  end
  iszero(total) || return SimplexSlowOverlap(SimplexOutsideSlowSector, zero_coefficient)

  coefficient = one(imaginary)
  prefix = zero(H)
  for index in 1:(N - 1)
    prefix += difference[index]
    iszero(prefix) && return SimplexSlowOverlap(SimplexSlowReducible, zero_coefficient)
    coefficient *= imaginary * inverse_weight(prefix)
  end
  return SimplexSlowOverlap(SimplexSlowPrimitive, coefficient)
end

function simplex_output_word_count(channel_count::Int, harmonic_count::Int, output_number::Int)
  channel_count >= 1 || throw(ArgumentError("channel count must be positive"))
  harmonic_count >= 1 || throw(ArgumentError("harmonic count must be positive"))
  output_number >= 0 || throw(ArgumentError("output number must be nonnegative"))
  return (channel_count * harmonic_count)^output_number
end

function simplex_output_triangle_count(channel_count::Int, harmonic_count::Int, order::Int)
  order >= 0 || throw(ArgumentError("order must be nonnegative"))
  return sum(
    simplex_output_word_count(channel_count, harmonic_count, output_number) for
    output_number in 0:(order + 1)
  )
end

function simplex_two_event_slow_paste(
  amplitudes::AbstractDict{H,T};
  channel,
  product,
  paste,
  imaginary,
  inverse_weight,
  zero_component,
  mismatch_abs=nothing,
) where {H,T}
  result = zero_component
  channels = (channel, channel)
  for (first_left, first_left_value) in amplitudes,
    (second_left, second_left_value) in amplitudes,
    (first_right, first_right_value) in amplitudes,
    (second_right, second_right_value) in amplitudes

    mismatch = first_left - first_right
    isnothing(mismatch_abs) || abs(mismatch) == mismatch_abs || continue

    left = SimplexOutputWord(channels, (first_left, second_left))
    right = SimplexOutputWord(channels, (first_right, second_right))
    overlap = simplex_slow_overlap(
      left, right; imaginary, inverse_weight
    )
    overlap.class == SimplexSlowPrimitive || continue

    left_value = product(second_left_value, first_left_value)
    right_value = product(second_right_value, first_right_value)
    result += overlap.coefficient * paste(left_value, right_value)
  end
  return result
end

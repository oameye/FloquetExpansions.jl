struct WaveOperator{H,P}
  coefficients::Vector{Dict{H,P}}
  bloch_effective_generator::Vector{P}
end

function series_entry!(series::Dict{H,P}, harmonic::H) where {H,P}
  haskey(series, harmonic) || (series[harmonic] = zero(P))
  return series[harmonic]
end

function add_periodic_scaled!(
  target::Dict{H,P}, source::Dict{H,P}, scale::Number
) where {H,P}
  for (harmonic, polynomial) in source
    add_scaled!(series_entry!(target, harmonic), polynomial, scale)
  end
  return target
end

function add_periodic_product!(
  target::Dict{H,P}, left::Dict{H,P}, right::Dict{H,P}, scale::Number
) where {H,P}
  for (left_harmonic, left_polynomial) in left, (right_harmonic, right_polynomial) in right
    entry = series_entry!(target, left_harmonic + right_harmonic)
    add_product!(entry, left_polynomial, right_polynomial, scale)
  end
  return target
end

function add_right_static_product!(
  target::Dict{H,P}, series::Dict{H,P}, static::P, scale::Number
) where {H,P}
  iszero(static) && return target
  for (harmonic, polynomial) in series
    add_product!(series_entry!(target, harmonic), polynomial, static, scale)
  end
  return target
end

bloch_harmonic_inverse(harmonic::Int) = 1 // harmonic

function bloch_antiderivative(residual::Dict{H,P}, zero_harmonic::H) where {H,P}
  wave = Dict{H,P}()
  for (harmonic, polynomial) in residual
    (harmonic == zero_harmonic || iszero(polynomial)) && continue
    wave[harmonic] = add_scaled!(zero(P), polynomial, bloch_harmonic_inverse(harmonic))
  end
  return wave
end

function bloch_residual(
  letters::Dict{H,P}, wave::Vector{Dict{H,P}}, bloch_effective_generator::Vector{P}, n::Int
) where {H,P}
  residual = add_periodic_product!(Dict{H,P}(), letters, wave[n], 1)
  for j in 1:n
    effective = bloch_effective_generator[n - j + 1]
    add_right_static_product!(residual, wave[j], effective, -1)
  end
  return residual
end

function bloch_wave_operator(letters::Dict{H,P}, order::Int, zero_harmonic::H) where {H,P}
  bloch_effective_generator = P[get(letters, zero_harmonic, zero(P))]
  wave = Dict{H,P}[]
  order > 1 && push!(wave, bloch_antiderivative(letters, zero_harmonic))
  for n in 1:(order - 1)
    residual = bloch_residual(letters, wave, bloch_effective_generator, n)
    push!(bloch_effective_generator, get(residual, zero_harmonic, zero(P)))
    n < order - 1 && push!(wave, bloch_antiderivative(residual, zero_harmonic))
  end
  return WaveOperator(wave, bloch_effective_generator)
end

function wave_times_static_factor(
  wave::Vector{Dict{H,P}}, static_factor::Vector{P}, n::Int
) where {H,P}
  prefactor = add_periodic_scaled!(Dict{H,P}(), wave[n], 1)
  for j in 1:(n - 1)
    add_right_static_product!(prefactor, wave[j], static_factor[n - j + 1], 1)
  end
  return prefactor
end

function log_series_nonlinear_terms!(
  powers::Dict{Tuple{Int,Int},Dict{H,P}}, normalized::Vector{Dict{H,P}}, n::Int
) where {H,P}
  nonlinear = Dict{H,P}()
  for power in 2:n
    power_coefficient = Dict{H,P}()
    for k in 1:(n - power + 1)
      previous = powers[(power - 1, n - k)]
      add_periodic_product!(power_coefficient, normalized[k], previous, 1)
    end
    powers[(power, n)] = power_coefficient
    add_periodic_scaled!(nonlinear, power_coefficient, (-1)^(power + 1) // power)
  end
  return nonlinear
end

function connected_log_series(
  wave::Vector{Dict{H,P}}, highest_power::Int, zero_harmonic::H
) where {H,P}
  static_factor = P[one(P)]
  normalized = Dict{H,P}[]
  powers = Dict{Tuple{Int,Int},Dict{H,P}}()
  connected_log = Dict{H,P}[]

  for n in 1:highest_power
    normalized_n = wave_times_static_factor(wave, static_factor, n)
    nonlinear = log_series_nonlinear_terms!(powers, normalized, n)
    static_part = get(nonlinear, zero_harmonic, zero(P))
    push!(static_factor, add_scaled!(zero(P), static_part, -1))
    add_scaled!(series_entry!(normalized_n, zero_harmonic), static_part, -1)
    push!(normalized, normalized_n)
    powers[(1, n)] = normalized_n
    connected_log_n = add_periodic_scaled!(Dict{H,P}(), normalized_n, 1)
    push!(connected_log, add_periodic_scaled!(connected_log_n, nonlinear, 1))
  end

  return (; connected_log, static_factor)
end

function static_series_inverse(series::Vector{P}, highest_power::Int) where {P}
  inverse = P[one(P)]
  for n in 1:highest_power
    coefficient = zero(P)
    for k in 1:n
      add_product!(coefficient, series[k + 1], inverse[n - k + 1], -1)
    end
    push!(inverse, coefficient)
  end
  return inverse
end

function static_series_product(
  left::Vector{P}, right::Vector{P}, highest_power::Int
) where {P}
  result = P[]
  for n in 0:highest_power
    coefficient = zero(P)
    for k in 0:n
      add_product!(coefficient, left[k + 1], right[n - k + 1], 1)
    end
    push!(result, coefficient)
  end
  return result
end

function normalize_to_van_vleck(wave::WaveOperator{H,P}, zero_harmonic::H) where {H,P}
  bloch_effective = wave.bloch_effective_generator
  highest_power = length(bloch_effective) - 1
  (; connected_log, static_factor) = connected_log_series(
    wave.coefficients, highest_power, zero_harmonic
  )
  inverse_static_factor = static_series_inverse(static_factor, highest_power)
  bloch_effective_times_static_factor = static_series_product(
    bloch_effective, static_factor, highest_power
  )
  effective = static_series_product(
    inverse_static_factor, bloch_effective_times_static_factor, highest_power
  )
  return (; connected_log, effective)
end

function van_vleck_words(letters::Vector{H}, order::Int, zero_harmonic::H) where {H}
  C = typeof(bloch_harmonic_inverse(one(H)))
  P = WordPolynomial{H,C}
  components = Dict{H,P}(letter => word_letter(letter, C) for letter in letters)
  wave = bloch_wave_operator(components, order, zero_harmonic)
  return normalize_to_van_vleck(wave, zero_harmonic)
end

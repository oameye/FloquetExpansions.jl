struct LyndonBracketNode
  left::Int
  right::Int
end

struct LyndonTerm{C}
  node::Int
  coefficient::C
end

struct ConnectedLogPlan{H,C}
  letters::Vector{H}
  brackets::Vector{LyndonBracketNode}
  outputs::Vector{Dict{H,Vector{LyndonTerm{C}}}}
end

function periodic_sum(left::Dict{H,T}, right::Dict{H,T}) where {H,T}
  result = copy(left)
  for (harmonic, value) in right
    result[harmonic] = haskey(result, harmonic) ? result[harmonic] + value : value
  end
  return result
end

function periodic_scale(scale::Number, series::Dict{H,T}) where {H,T}
  return Dict{H,T}(harmonic => scale * value for (harmonic, value) in series)
end

function periodic_product(left::Dict{H,T}, right::Dict{H,T}) where {H,T}
  result = Dict{H,T}()
  for (left_harmonic, left_value) in left, (right_harmonic, right_value) in right
    harmonic = left_harmonic + right_harmonic
    value = left_value * right_value
    result[harmonic] = haskey(result, harmonic) ? result[harmonic] + value : value
  end
  return result
end

function right_static_product(series::Dict{H,T}, static::T) where {H,T}
  return Dict{H,T}(harmonic => value * static for (harmonic, value) in series)
end

function wave_times_static_factor(
  wave::Vector{Dict{H,T}}, static_factor::Vector{T}, n::Int
) where {H,T}
  result = copy(wave[n])
  for j in 1:(n - 1)
    result = periodic_sum(result, right_static_product(wave[j], static_factor[n - j + 1]))
  end
  return result
end

function log_series_nonlinear_terms!(
  powers::Dict{Tuple{Int,Int},Dict{H,T}}, normalized::Vector{Dict{H,T}}, n::Int
) where {H,T}
  nonlinear = Dict{H,T}()
  for power in 2:n
    power_coefficient = Dict{H,T}()
    for k in 1:(n - power + 1)
      term = periodic_product(normalized[k], powers[(power - 1, n - k)])
      power_coefficient = periodic_sum(power_coefficient, term)
    end
    powers[(power, n)] = power_coefficient
    log_weight = (-1)^(power + 1) // power
    nonlinear = periodic_sum(nonlinear, periodic_scale(log_weight, power_coefficient))
  end
  return nonlinear
end

function connected_log_series(
  wave::Vector{Dict{H,T}}, highest_power::Int, zero_harmonic::H
) where {H,T}
  static_factor = T[one(T)]
  normalized = Dict{H,T}[]
  powers = Dict{Tuple{Int,Int},Dict{H,T}}()
  connected_log = Dict{H,T}[]

  for n in 1:highest_power
    prefactor = wave_times_static_factor(wave, static_factor, n)
    nonlinear = log_series_nonlinear_terms!(powers, normalized, n)
    static_part = get(periodic_sum(prefactor, nonlinear), zero_harmonic, zero(one(T)))
    push!(static_factor, -static_part)

    normalized_n = periodic_sum(prefactor, Dict(zero_harmonic => -static_part))
    push!(normalized, normalized_n)
    powers[(1, n)] = normalized_n
    push!(connected_log, periodic_sum(normalized_n, nonlinear))
  end

  return connected_log
end

function connected_log_words(plan::WaveOperatorPlan{H}) where {H}
  C = typeof(bloch_harmonic_inverse(one(H)))
  P = WordPolynomial{H,C}
  components = Dict{H,P}(letter => word_letter(letter, C) for letter in plan.support)
  wave = evaluate_wave_operator(plan, components, P(), BlochConventions(*, 1, 1))
  return connected_log_series(wave.coefficients, plan.order - 1, plan.zero_harmonic)
end

function lyndon_node!(
  word::Vector{H},
  letter_ids::Dict{H,Int},
  word_ids::Dict{Vector{H},Int},
  brackets::Vector{LyndonBracketNode},
) where {H}
  length(word) == 1 && return letter_ids[only(word)]
  haskey(word_ids, word) && return word_ids[word]
  prefix, suffix = lyndon_standard_factorization(word)
  left = lyndon_node!(prefix, letter_ids, word_ids, brackets)
  right = lyndon_node!(suffix, letter_ids, word_ids, brackets)
  push!(brackets, LyndonBracketNode(left, right))
  node = length(letter_ids) + length(brackets)
  word_ids[word] = node
  return node
end

function compile_connected_log_plan(plan::WaveOperatorPlan{H}) where {H}
  letters = plan.support
  log_words = connected_log_words(plan)
  C = typeof(bloch_harmonic_inverse(one(H)))
  letter_ids = Dict(letter => index for (index, letter) in enumerate(letters))
  word_ids = Dict{Vector{H},Int}()
  brackets = LyndonBracketNode[]
  cache = Dict{Vector{H},WordPolynomial{H,C}}()
  outputs = Dict{H,Vector{LyndonTerm{C}}}[]

  for series in log_words
    compiled = Dict{H,Vector{LyndonTerm{C}}}()
    for (harmonic, polynomial) in series
      iszero(polynomial) && continue
      compiled[harmonic] = LyndonTerm{C}[
        LyndonTerm(lyndon_node!(word, letter_ids, word_ids, brackets), coefficient) for
        (word, coefficient) in lyndon_coordinates(polynomial, cache)
      ]
    end
    push!(outputs, compiled)
  end

  return ConnectedLogPlan(letters, brackets, outputs)
end

function evaluate_connected_log(
  plan::ConnectedLogPlan{H,C}, components::AbstractDict{H,T}, zero_component::T, product
) where {H,C,T}
  nletters = length(plan.letters)
  values = Vector{T}(undef, nletters + length(plan.brackets))
  for (index, letter) in enumerate(plan.letters)
    values[index] = components[letter]
  end
  for (index, bracket) in enumerate(plan.brackets)
    left = values[bracket.left]
    right = values[bracket.right]
    values[nletters + index] =
      simplify_component(product(left, right) - product(right, left))::T
  end

  connected_log = Dict{H,T}[]
  for output in plan.outputs
    series = Dict{H,T}()
    for (harmonic, terms) in output
      value = zero_component
      for term in terms
        value += term.coefficient * values[term.node]
      end
      component = simplify_component(value)::T
      iszero(component) || (series[harmonic] = component)
    end
    push!(connected_log, series)
  end
  return connected_log
end

function apply_weight_phase!(connected_log::Vector{Dict{H,T}}, weight_phase) where {H,T}
  for (n, series) in enumerate(connected_log)
    scale = weight_phase^n
    for harmonic in collect(keys(series))
      component = simplify_component(scale * series[harmonic])::T
      iszero(component) ? delete!(series, harmonic) : (series[harmonic] = component)
    end
  end
  return connected_log
end

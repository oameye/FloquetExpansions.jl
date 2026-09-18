function sparse_accumulate!(terms::Dict{Tuple,C}, word::Tuple, coefficient) where {C}
  value = get(terms, word, zero(C)) + convert(C, coefficient)
  if iszero(value)
    haskey(terms, word) && delete!(terms, word)
  else
    terms[word] = value
  end
  return terms
end

function sparse_word_terms!(embedding::Dict{H,Dict{Tuple,C}}, harmonic::H) where {H,C}
  return get!(embedding, harmonic) do
    return Dict{Tuple,C}()
  end
end

function sparse_add_terms!(
  destination::Dict{Tuple,C}, source::Dict{Tuple,C}, weight::C=one(C)
) where {C}
  for (word, coefficient) in source
    sparse_accumulate!(destination, word, weight * coefficient)
  end
  return destination
end

function sparse_copy_embedding(embedding::Dict{H,Dict{Tuple,C}}) where {H,C}
  return Dict(harmonic => copy(terms) for (harmonic, terms) in embedding)
end

function sparse_add_embedding!(
  destination::Dict{H,Dict{Tuple,C}}, source::Dict{H,Dict{Tuple,C}}, weight::C=one(C)
) where {H,C}
  for (harmonic, terms) in source
    sparse_add_terms!(sparse_word_terms!(destination, harmonic), terms, weight)
  end
  return destination
end

function sparse_periodic_product(
  left::Dict{H,Dict{Tuple,C}}, right::Dict{H,Dict{Tuple,C}}
) where {H,C}
  result = Dict{H,Dict{Tuple,C}}()
  for (left_harmonic, left_terms) in left, (right_harmonic, right_terms) in right
    output = sparse_word_terms!(result, left_harmonic + right_harmonic)
    for (left_word, left_coefficient) in left_terms,
      (right_word, right_coefficient) in right_terms

      sparse_accumulate!(
        output, (left_word..., right_word...), left_coefficient * right_coefficient
      )
    end
  end
  return result
end

function sparse_right_static_product(
  periodic::Dict{H,Dict{Tuple,C}}, static::Dict{Tuple,C}
) where {H,C}
  result = Dict{H,Dict{Tuple,C}}()
  for (harmonic, periodic_terms) in periodic
    output = sparse_word_terms!(result, harmonic)
    for (periodic_word, periodic_coefficient) in periodic_terms,
      (static_word, static_coefficient) in static

      sparse_accumulate!(
        output,
        (periodic_word..., static_word...),
        periodic_coefficient * static_coefficient,
      )
    end
  end
  return result
end

function sparse_scale_embedding!(embedding::Dict{H,Dict{Tuple,C}}, weight::C) where {H,C}
  for terms in values(embedding)
    for word in collect(keys(terms))
      coefficient = weight * terms[word]
      if iszero(coefficient)
        delete!(terms, word)
      else
        terms[word] = coefficient
      end
    end
  end
  return embedding
end

function compile_sparse_connected_log_words(support_input, order::Int)
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  support = sort(unique(collect(support_input)))
  isempty(support) && throw(ArgumentError("harmonic support must not be empty"))

  H = eltype(support)
  C = Rational{Int}
  zero_harmonic = zero(first(support))
  one_coefficient = one(C)

  effective = Vector{Dict{Tuple,C}}()
  B0 = Dict{Tuple,C}()
  if zero_harmonic in support
    B0[(zero_harmonic,)] = one_coefficient
  end
  push!(effective, B0)

  wave = Vector{Dict{H,Dict{Tuple,C}}}()
  if order > 1
    X1 = Dict{H,Dict{Tuple,C}}()
    for harmonic in support
      iszero(harmonic) && continue
      sparse_word_terms!(X1, harmonic)[(harmonic,)] = 1 // harmonic
    end
    push!(wave, X1)
  end

  for n in 1:(order - 1)
    residual = Dict{H,Dict{Tuple,C}}()

    for generator_harmonic in support, (wave_harmonic, wave_terms) in wave[n]
      output = sparse_word_terms!(residual, generator_harmonic + wave_harmonic)
      for (wave_word, wave_coefficient) in wave_terms
        sparse_accumulate!(output, (generator_harmonic, wave_word...), wave_coefficient)
      end
    end

    for wave_order in 1:n
      effective_order = n - wave_order
      B = effective[effective_order + 1]
      isempty(B) && continue
      for (harmonic, wave_terms) in wave[wave_order]
        output = sparse_word_terms!(residual, harmonic)
        for (wave_word, wave_coefficient) in wave_terms,
          (effective_word, effective_coefficient) in B

          sparse_accumulate!(
            output,
            (wave_word..., effective_word...),
            -wave_coefficient * effective_coefficient,
          )
        end
      end
    end

    push!(effective, copy(get(residual, zero_harmonic, Dict{Tuple,C}())))

    if n < order - 1
      Xnext = Dict{H,Dict{Tuple,C}}()
      for (harmonic, residual_terms) in residual
        iszero(harmonic) && continue
        output = sparse_word_terms!(Xnext, harmonic)
        weight = 1 // harmonic
        for (word, coefficient) in residual_terms
          sparse_accumulate!(output, word, weight * coefficient)
        end
      end
      push!(wave, Xnext)
    end
  end

  canonical_order = order - 1
  static_factor = Dict{Tuple,C}[Dict(() => one_coefficient)]
  normalized = Vector{Dict{H,Dict{Tuple,C}}}()
  log_embedding = Vector{Dict{H,Dict{Tuple,C}}}()
  powers = [
    [Dict{H,Dict{Tuple,C}}() for _ in 1:max(canonical_order, 1)] for
    _ in 1:max(canonical_order, 1)
  ]

  for n in 1:canonical_order
    prefactor = sparse_copy_embedding(wave[n])
    for j in 1:(n - 1)
      correction = sparse_right_static_product(wave[j], static_factor[n - j + 1])
      sparse_add_embedding!(prefactor, correction)
    end

    nonlinear_log = Dict{H,Dict{Tuple,C}}()
    for power in 2:n
      power_coefficient = Dict{H,Dict{Tuple,C}}()
      for k in 1:(n - power + 1)
        contribution = sparse_periodic_product(normalized[k], powers[power - 1][n - k])
        sparse_add_embedding!(power_coefficient, contribution)
      end
      powers[power][n] = power_coefficient
      weight = convert(C, (-1)^(power + 1) * (1 // power))
      contribution = sparse_copy_embedding(power_coefficient)
      sparse_scale_embedding!(contribution, weight)
      sparse_add_embedding!(nonlinear_log, contribution)
    end

    candidate = sparse_copy_embedding(prefactor)
    sparse_add_embedding!(candidate, nonlinear_log)
    static_n = Dict{Tuple,C}()
    for (word, coefficient) in get(candidate, zero_harmonic, Dict{Tuple,C}())
      sparse_accumulate!(static_n, word, -coefficient)
    end
    push!(static_factor, static_n)

    normalized_n = sparse_copy_embedding(prefactor)
    sparse_add_terms!(sparse_word_terms!(normalized_n, zero_harmonic), static_n)
    push!(normalized, normalized_n)
    powers[1][n] = normalized_n

    log_n = sparse_copy_embedding(normalized_n)
    sparse_add_embedding!(log_n, nonlinear_log)
    push!(log_embedding, log_n)
  end

  return log_embedding
end

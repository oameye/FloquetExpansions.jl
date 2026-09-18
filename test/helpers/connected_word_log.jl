using FloquetExpansions

const FE_CWL = FloquetExpansions

struct HarmonicWordPolynomial{C}
  terms::Dict{Tuple,C}
end

HarmonicWordPolynomial{C}() where {C} = HarmonicWordPolynomial{C}(Dict{Tuple,C}())

function harmonic_word_leaf(harmonic, ::Type{C}=Rational{Int}) where {C}
  return HarmonicWordPolynomial{C}(Dict((harmonic,) => one(C)))
end

function harmonic_word_polynomial(terms::Dict{Tuple,C}) where {C}
  cleaned = Dict{Tuple,C}()
  for (word, coefficient) in terms
    iszero(coefficient) || (cleaned[word] = coefficient)
  end
  return HarmonicWordPolynomial{C}(cleaned)
end

Base.zero(polynomial::HarmonicWordPolynomial{C}) where {C} = HarmonicWordPolynomial{C}()
function Base.one(polynomial::HarmonicWordPolynomial{C}) where {C}
  return HarmonicWordPolynomial{C}(Dict(() => one(C)))
end
Base.iszero(polynomial::HarmonicWordPolynomial) = isempty(polynomial.terms)

function Base.:+(
  left::HarmonicWordPolynomial{C}, right::HarmonicWordPolynomial{C}
) where {C}
  terms = copy(left.terms)
  for (word, coefficient) in right.terms
    terms[word] = get(terms, word, zero(C)) + coefficient
  end
  return harmonic_word_polynomial(terms)
end

function Base.:-(polynomial::HarmonicWordPolynomial{C}) where {C}
  return HarmonicWordPolynomial{C}(
    Dict(word => -coefficient for (word, coefficient) in polynomial.terms)
  )
end
function Base.:-(
  left::HarmonicWordPolynomial{C}, right::HarmonicWordPolynomial{C}
) where {C}
  return left + (-right)
end

function Base.:*(weight::Number, polynomial::HarmonicWordPolynomial{C}) where {C}
  terms = Dict{Tuple,C}()
  for (word, coefficient) in polynomial.terms
    scaled = convert(C, weight * coefficient)
    iszero(scaled) || (terms[word] = scaled)
  end
  return HarmonicWordPolynomial{C}(terms)
end

Base.:*(polynomial::HarmonicWordPolynomial, weight::Number) = weight * polynomial

function harmonic_word_product(
  left::HarmonicWordPolynomial{C}, right::HarmonicWordPolynomial{C}
) where {C}
  terms = Dict{Tuple,C}()
  for (left_word, left_coefficient) in left.terms,
    (right_word, right_coefficient) in right.terms

    word = (left_word..., right_word...)
    terms[word] = get(terms, word, zero(C)) + left_coefficient * right_coefficient
  end
  return harmonic_word_polynomial(terms)
end

function harmonic_word_commutator(
  left::HarmonicWordPolynomial{C}, right::HarmonicWordPolynomial{C}
) where {C}
  return harmonic_word_product(left, right) - harmonic_word_product(right, left)
end

function dynkin_word(word::Tuple, ::Type{C}) where {C}
  isempty(word) && return HarmonicWordPolynomial{C}()
  result = harmonic_word_leaf(first(word), C)
  for harmonic in Iterators.drop(word, 1)
    result = harmonic_word_commutator(result, harmonic_word_leaf(harmonic, C))
  end
  return result
end

function dynkin_projection(polynomial::HarmonicWordPolynomial{C}) where {C}
  result = HarmonicWordPolynomial{C}()
  for (word, coefficient) in polynomial.terms
    result += coefficient * dynkin_word(word, C)
  end
  return result
end

function harmonic_word_lengths(polynomial::HarmonicWordPolynomial)
  return unique(length(word) for word in keys(polynomial.terms))
end

function is_lyndon_word(word::Tuple)
  isempty(word) && return false
  return all(isless(word, word[index:end]) for index in 2:length(word))
end

function lyndon_standard_factorization(word::Tuple)
  length(word) > 1 || throw(ArgumentError("a Lyndon leaf has no standard factorization"))
  for split in 2:length(word)
    suffix = word[split:end]
    if is_lyndon_word(suffix)
      return word[1:(split - 1)], suffix
    end
  end
  return throw(ArgumentError("word is not Lyndon"))
end

function lyndon_bracket_polynomial!(
  cache::Dict{Tuple,HarmonicWordPolynomial{C}}, word::Tuple, ::Type{C}
) where {C}
  haskey(cache, word) && return cache[word]
  result = if length(word) == 1
    harmonic_word_leaf(first(word), C)
  else
    left_word, right_word = lyndon_standard_factorization(word)
    left = lyndon_bracket_polynomial!(cache, left_word, C)
    right = lyndon_bracket_polynomial!(cache, right_word, C)
    harmonic_word_commutator(left, right)
  end
  cache[word] = result
  return result
end

function lyndon_bracket_polynomial(word::Tuple, ::Type{C}) where {C}
  cache = Dict{Tuple,HarmonicWordPolynomial{C}}()
  return lyndon_bracket_polynomial!(cache, word, C)
end

function lyndon_decomposition(
  polynomial::HarmonicWordPolynomial{C},
  bracket_cache::Dict{Tuple,HarmonicWordPolynomial{C}},
) where {C}
  residual = copy(polynomial.terms)
  coefficients = Dict{Tuple,C}()
  while !isempty(residual)
    word = minimum(keys(residual))
    is_lyndon_word(word) ||
      throw(ArgumentError("primitive polynomial has a non-Lyndon leading word"))
    coefficient = residual[word]
    coefficients[word] = get(coefficients, word, zero(C)) + coefficient
    bracket = lyndon_bracket_polynomial!(bracket_cache, word, C)
    for (term_word, term_coefficient) in bracket.terms
      updated = get(residual, term_word, zero(C)) - coefficient * term_coefficient
      if iszero(updated)
        delete!(residual, term_word)
      else
        residual[term_word] = updated
      end
    end
  end
  return coefficients
end

function lyndon_decomposition(polynomial::HarmonicWordPolynomial{C}) where {C}
  cache = Dict{Tuple,HarmonicWordPolynomial{C}}()
  return lyndon_decomposition(polynomial, cache)
end

function collect_lyndon_bracket_nodes!(nodes::Set{Tuple}, word::Tuple)
  length(word) <= 1 && return nodes
  push!(nodes, word)
  left, right = lyndon_standard_factorization(word)
  collect_lyndon_bracket_nodes!(nodes, left)
  collect_lyndon_bracket_nodes!(nodes, right)
  return nodes
end

function lyndon_profile(embeddings::AbstractVector)
  coefficient_count = 0
  basis_words = Set{Tuple}()
  bracket_nodes = Set{Tuple}()
  bracket_cache = Dict{Tuple,HarmonicWordPolynomial{Rational{Int}}}()
  for embedding in embeddings, polynomial in values(embedding)
    decomposition = lyndon_decomposition(polynomial, bracket_cache)
    coefficient_count += length(decomposition)
    for word in keys(decomposition)
      push!(basis_words, word)
      collect_lyndon_bracket_nodes!(bracket_nodes, word)
    end
  end
  return (
    coefficients=coefficient_count,
    basis_words=length(basis_words),
    bracket_nodes=length(bracket_nodes),
    associative_products=2 * length(bracket_nodes),
  )
end

function evaluate_harmonic_word_polynomial(
  polynomial::HarmonicWordPolynomial,
  components::AbstractDict;
  product,
  zero_component,
  simplifier=identity,
)
  isempty(polynomial.terms) && return zero_component
  identity_component = one(first(values(components)))
  result = zero_component
  for (word, coefficient) in polynomial.terms
    value = identity_component
    for harmonic in word
      value = product(value, components[harmonic])
    end
    result += coefficient * value
  end
  return simplifier(result)
end

function evaluate_harmonic_word_embedding(
  embedding::Dict{H,<:HarmonicWordPolynomial},
  components::AbstractDict{H,T};
  product,
  zero_component::T,
  simplifier=identity,
) where {H,T}
  result = Dict{H,T}()
  for (harmonic, polynomial) in embedding
    value = evaluate_harmonic_word_polynomial(
      polynomial, components; product, zero_component, simplifier
    )
    iszero(value) || (result[harmonic] = value)
  end
  return result
end

function connected_word_reconstruction(support, order::Int)
  C = Rational{Int}
  zero_harmonic = zero(first(support))
  components = Dict(harmonic => harmonic_word_leaf(harmonic, C) for harmonic in support)
  zero_component = HarmonicWordPolynomial{C}()
  plan = FE_CWL.compile_bloch_projection_plan(support, order, zero_harmonic)
  bloch = FE_CWL.evaluate_bloch_projection_plan(
    plan,
    components;
    product=harmonic_word_product,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component,
  )
  converted = FE_CWL.bloch_van_vleck_reconstruction(
    plan, bloch; product=harmonic_word_product, zero_component
  )
  return plan, bloch, converted
end

function connected_word_log_only(support, order::Int)
  C = Rational{Int}
  zero_harmonic = zero(first(support))
  components = Dict(harmonic => harmonic_word_leaf(harmonic, C) for harmonic in support)
  zero_component = HarmonicWordPolynomial{C}()
  plan = FE_CWL.compile_bloch_projection_plan(support, order, zero_harmonic)
  bloch = FE_CWL.evaluate_bloch_projection_plan(
    plan,
    components;
    product=harmonic_word_product,
    inverse_weight=harmonic -> 1 // harmonic,
    zero_component,
  )

  reconstruction_order = length(bloch.effective) - 1
  identity_component = one(first(bloch.effective))
  static_factor = typeof(identity_component)[identity_component]
  normalized_embedding = Vector{Dict{typeof(zero_harmonic),HarmonicWordPolynomial{C}}}()
  log_embedding = Vector{Dict{typeof(zero_harmonic),HarmonicWordPolynomial{C}}}()
  powers = [
    [Dict{typeof(zero_harmonic),HarmonicWordPolynomial{C}}() for _ in 1:max(reconstruction_order, 1)] for _ in 1:max(reconstruction_order, 1)
  ]
  counts = FE_CWL.BlochVanVleckCounts()

  for n in 1:reconstruction_order
    prefactor = copy(bloch.wave[n])
    for j in 1:(n - 1)
      counts.factor_products += 1
      correction = FE_CWL.bloch_vv_right_static_product(
        bloch.wave[j], static_factor[n - j + 1], harmonic_word_product, identity, counts
      )
      prefactor = FE_CWL.bloch_vv_add(prefactor, correction, identity)
    end

    nonlinear_log = Dict{typeof(zero_harmonic),HarmonicWordPolynomial{C}}()
    for power in 2:n
      power_coefficient = Dict{typeof(zero_harmonic),HarmonicWordPolynomial{C}}()
      for k in 1:(n - power + 1)
        counts.log_products += 1
        contribution = FE_CWL.bloch_vv_periodic_product(
          normalized_embedding[k],
          powers[power - 1][n - k],
          harmonic_word_product,
          identity,
          counts,
        )
        power_coefficient = FE_CWL.bloch_vv_add(power_coefficient, contribution, identity)
      end
      powers[power][n] = power_coefficient
      weight = (-1)^(power + 1) * (1 // power)
      nonlinear_log = FE_CWL.bloch_vv_add(
        nonlinear_log, FE_CWL.bloch_vv_scale(weight, power_coefficient, identity), identity
      )
    end

    candidate = FE_CWL.bloch_vv_add(prefactor, nonlinear_log, identity)
    static_n = -get(candidate, zero_harmonic, zero_component)
    push!(static_factor, static_n)

    normalized_n = copy(prefactor)
    if !iszero(static_n)
      FE_CWL.bloch_vv_accumulate!(normalized_n, zero_harmonic, static_n)
      normalized_n = FE_CWL.bloch_vv_simplify_embedding(normalized_n, identity)
    end
    push!(normalized_embedding, normalized_n)
    powers[1][n] = normalized_n
    push!(log_embedding, FE_CWL.bloch_vv_add(normalized_n, nonlinear_log, identity))
  end

  return log_embedding
end

function harmonic_word_count(embedding)
  return sum(length(polynomial.terms) for polynomial in values(embedding))
end

function harmonic_word_count(embeddings::AbstractVector)
  return sum(harmonic_word_count, embeddings; init=0)
end

function harmonic_word_prefix_products(embeddings::AbstractVector)
  prefixes = Set{Tuple}()
  for embedding in embeddings,
    polynomial in values(embedding),
    word in keys(polynomial.terms)

    for length_prefix in 2:length(word)
      push!(prefixes, word[1:length_prefix])
    end
  end
  return length(prefixes)
end

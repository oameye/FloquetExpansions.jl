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

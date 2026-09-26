struct WordPolynomial{H,C}
  terms::Dict{Vector{H},C}
end

WordPolynomial{H,C}() where {H,C} = WordPolynomial(Dict{Vector{H},C}())

function word_letter(letter::H, ::Type{C}) where {H,C}
  return WordPolynomial(Dict{Vector{H},C}(H[letter] => one(C)))
end

Base.one(::Type{WordPolynomial{H,C}}) where {H,C} = WordPolynomial(Dict(H[] => one(C)))
Base.zero(::WordPolynomial{H,C}) where {H,C} = WordPolynomial{H,C}()
Base.iszero(polynomial::WordPolynomial) = isempty(polynomial.terms)

simplify_component(polynomial::WordPolynomial) = polynomial

function accumulate_word!(
  terms::Dict{Vector{H},C}, word::Vector{H}, coefficient
) where {H,C}
  value = get(terms, word, zero(C)) + convert(C, coefficient)
  iszero(value) ? delete!(terms, word) : (terms[word] = value)
  return terms
end

function Base.:+(left::WordPolynomial{H,C}, right::WordPolynomial{H,C}) where {H,C}
  terms = copy(left.terms)
  for (word, coefficient) in right.terms
    accumulate_word!(terms, word, coefficient)
  end
  return WordPolynomial(terms)
end

function Base.:-(left::WordPolynomial{H,C}, right::WordPolynomial{H,C}) where {H,C}
  terms = copy(left.terms)
  for (word, coefficient) in right.terms
    accumulate_word!(terms, word, -coefficient)
  end
  return WordPolynomial(terms)
end

Base.:-(polynomial::WordPolynomial) = -1 * polynomial

function Base.:*(scale::Number, polynomial::WordPolynomial{H,C}) where {H,C}
  terms = Dict{Vector{H},C}()
  for (word, coefficient) in polynomial.terms
    accumulate_word!(terms, word, scale * coefficient)
  end
  return WordPolynomial(terms)
end

function Base.:*(left::WordPolynomial{H,C}, right::WordPolynomial{H,C}) where {H,C}
  terms = Dict{Vector{H},C}()
  for (left_word, left_coefficient) in left.terms,
    (right_word, right_coefficient) in right.terms

    accumulate_word!(
      terms, vcat(left_word, right_word), left_coefficient * right_coefficient
    )
  end
  return WordPolynomial(terms)
end

word_commutator(left::WordPolynomial, right::WordPolynomial) = left * right - right * left

function is_lyndon_word(word::Vector)
  isempty(word) && return false
  return all(word < word[index:end] for index in 2:length(word))
end

function lyndon_standard_factorization(word::Vector)
  length(word) > 1 || throw(ArgumentError("a Lyndon letter has no standard factorization"))
  for split in 2:length(word)
    suffix = word[split:end]
    is_lyndon_word(suffix) && return word[1:(split - 1)], suffix
  end
  return throw(ArgumentError("word is not Lyndon"))
end

function lyndon_bracket!(
  cache::Dict{Vector{H},WordPolynomial{H,C}}, word::Vector{H}
) where {H,C}
  haskey(cache, word) && return cache[word]
  bracket = if length(word) == 1
    word_letter(only(word), C)
  else
    prefix, suffix = lyndon_standard_factorization(word)
    word_commutator(lyndon_bracket!(cache, prefix), lyndon_bracket!(cache, suffix))
  end
  cache[word] = bracket
  return bracket
end

function lyndon_coordinates(
  polynomial::WordPolynomial{H,C}, cache::Dict{Vector{H},WordPolynomial{H,C}}
) where {H,C}
  residual = copy(polynomial.terms)
  coordinates = Dict{Vector{H},C}()
  while !isempty(residual)
    word = minimum(keys(residual))
    is_lyndon_word(word) ||
      throw(ArgumentError("a Lie polynomial cannot have a non-Lyndon leading word"))
    coefficient = residual[word]
    coordinates[word] = coefficient
    for (term, term_coefficient) in lyndon_bracket!(cache, word).terms
      accumulate_word!(residual, term, -coefficient * term_coefficient)
    end
  end
  return coordinates
end

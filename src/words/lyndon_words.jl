struct WordPolynomial{H,C}
  terms::Dict{Vector{H},C}
end

WordPolynomial{H,C}() where {H,C} = WordPolynomial(Dict{Vector{H},C}())

function word_letter(letter::H, ::Type{C}) where {H,C}
  return WordPolynomial(Dict{Vector{H},C}(H[letter] => one(C)))
end

Base.one(::Type{WordPolynomial{H,C}}) where {H,C} = WordPolynomial(Dict(H[] => one(C)))
Base.zero(::Type{WordPolynomial{H,C}}) where {H,C} = WordPolynomial{H,C}()
Base.iszero(polynomial::WordPolynomial) = isempty(polynomial.terms)

function accumulate_word!(
  terms::Dict{Vector{H},C}, word::Vector{H}, coefficient
) where {H,C}
  value = get(terms, word, zero(C)) + convert(C, coefficient)
  iszero(value) ? delete!(terms, word) : (terms[word] = value)
  return terms
end

function add_scaled!(
  target::WordPolynomial{H,C}, source::WordPolynomial{H,C}, scale::Number
) where {H,C}
  for (word, coefficient) in source.terms
    accumulate_word!(target.terms, word, scale * coefficient)
  end
  return target
end

function add_product!(
  target::WordPolynomial{H,C},
  left::WordPolynomial{H,C},
  right::WordPolynomial{H,C},
  scale::Number,
) where {H,C}
  for (left_word, left_coefficient) in left.terms,
    (right_word, right_coefficient) in right.terms

    accumulate_word!(
      target.terms,
      vcat(left_word, right_word),
      scale * left_coefficient * right_coefficient,
    )
  end
  return target
end

function word_commutator(left::WordPolynomial{H,C}, right::WordPolynomial{H,C}) where {H,C}
  commutator = add_product!(zero(WordPolynomial{H,C}), left, right, 1)
  return add_product!(commutator, right, left, -1)
end

function is_lyndon_word(word::Vector)
  isempty(word) && return false
  return all(word < view(word, index:lastindex(word)) for index in 2:lastindex(word))
end

function lyndon_standard_factorization(word::Vector)
  length(word) > 1 || throw(ArgumentError("a Lyndon letter has no standard factorization"))
  for split in 2:(length(word) - 1)
    suffix = word[split:end]
    is_lyndon_word(suffix) && return word[1:(split - 1)], suffix
  end
  return word[1:(end - 1)], word[end:end]
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

function subtract_lyndon_bracket!(
  residual::Dict{Vector{H},C},
  cache::Dict{Vector{H},WordPolynomial{H,C}},
  word::Vector{H},
  coefficient::C,
) where {H,C}
  for (term, term_coefficient) in lyndon_bracket!(cache, word).terms
    accumulate_word!(residual, term, -coefficient * term_coefficient)
  end
  return residual
end

function lyndon_pass!(
  coordinates::Dict{Vector{H},C},
  residual::Dict{Vector{H},C},
  cache::Dict{Vector{H},WordPolynomial{H,C}},
) where {H,C}
  words = sort!(collect(keys(residual)))
  is_lyndon_word(first(words)) ||
    throw(ArgumentError("a Lie polynomial cannot have a non-Lyndon leading word"))
  for word in words
    (haskey(residual, word) && is_lyndon_word(word)) || continue
    coefficient = residual[word]
    accumulate_word!(coordinates, word, coefficient)
    subtract_lyndon_bracket!(residual, cache, word, coefficient)
  end
  return coordinates
end

function lyndon_coordinates(
  polynomial::WordPolynomial{H,C}, cache::Dict{Vector{H},WordPolynomial{H,C}}
) where {H,C}
  residual = copy(polynomial.terms)
  coordinates = Dict{Vector{H},C}()
  while !isempty(residual)
    lyndon_pass!(coordinates, residual, cache)
  end
  return coordinates
end

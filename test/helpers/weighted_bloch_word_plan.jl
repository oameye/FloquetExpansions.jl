mutable struct BlochWordPlanCounts
  generator_product_terms::Int
  folded_counterterms::Int
  coefficient_cancellations::Int
end

BlochWordPlanCounts() = BlochWordPlanCounts(0, 0, 0)

struct BlochWordPlan{H,C}
  support::Vector{H}
  order::Int
  wave::Vector{Dict{H,Dict{Tuple,C}}}
  effective::Vector{Dict{Tuple,C}}
  counts::BlochWordPlanCounts
end

function add_word_coefficient!(
  terms::Dict{Tuple,C}, word::Tuple, coefficient, counts::BlochWordPlanCounts
) where {C}
  value = get(terms, word, zero(C)) + convert(C, coefficient)
  if iszero(value)
    haskey(terms, word) && delete!(terms, word)
    counts.coefficient_cancellations += 1
  else
    terms[word] = value
  end
  return terms
end

function harmonic_word_terms!(
  terms::Dict{H,Dict{Tuple,C}}, harmonic::H
) where {H,C}
  return get!(terms, harmonic) do
    Dict{Tuple,C}()
  end
end

function compile_bloch_word_plan(
  support_input,
  order::Int;
  inverse_weight,
  zero_harmonic=zero(first(support_input)),
)
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  isempty(support_input) && throw(ArgumentError("harmonic support must not be empty"))

  support = unique(collect(support_input))
  H = typeof(zero_harmonic)
  all(harmonic isa H for harmonic in support) ||
    throw(ArgumentError("harmonic support and zero harmonic must have one common type"))
  support = H[harmonic for harmonic in support]

  nonzero_index = findfirst(!iszero, support)
  coefficient_sample =
    isnothing(nonzero_index) ? 1 : inverse_weight(support[nonzero_index])
  C = typeof(coefficient_sample)
  one_coefficient = one(coefficient_sample)

  wave = Vector{Dict{H,Dict{Tuple,C}}}()
  effective = Vector{Dict{Tuple,C}}()
  counts = BlochWordPlanCounts()

  B0 = Dict{Tuple,C}()
  if zero_harmonic in support
    B0[(zero_harmonic,)] = one_coefficient
  end
  push!(effective, B0)

  if order > 1
    X1 = Dict{H,Dict{Tuple,C}}()
    for harmonic in support
      iszero(harmonic) && continue
      harmonic_word_terms!(X1, harmonic)[(harmonic,)] =
        convert(C, inverse_weight(harmonic))
    end
    push!(wave, X1)
  end

  for n in 1:(order - 1)
    residual = Dict{H,Dict{Tuple,C}}()

    for (wave_harmonic, wave_terms) in wave[n], generator_harmonic in support
      output_harmonic = wave_harmonic + generator_harmonic
      output_terms = harmonic_word_terms!(residual, output_harmonic)
      for (wave_word, coefficient) in wave_terms
        counts.generator_product_terms += 1
        add_word_coefficient!(
          output_terms,
          (wave_word..., generator_harmonic),
          coefficient,
          counts,
        )
      end
    end

    for j in 1:n
      Bj = effective[n - j + 1]
      isempty(Bj) && continue
      for (wave_harmonic, wave_terms) in wave[j]
        output_terms = harmonic_word_terms!(residual, wave_harmonic)
        for (wave_word, wave_coefficient) in wave_terms,
          (effective_word, effective_coefficient) in Bj
          counts.folded_counterterms += 1
          add_word_coefficient!(
            output_terms,
            (effective_word..., wave_word...),
            -wave_coefficient * effective_coefficient,
            counts,
          )
        end
      end
    end

    Bn = get(residual, zero_harmonic, Dict{Tuple,C}())
    push!(effective, copy(Bn))

    if n < order - 1
      Xnext = Dict{H,Dict{Tuple,C}}()
      for (harmonic, residual_terms) in residual
        iszero(harmonic) && continue
        weight = inverse_weight(harmonic)
        output_terms = harmonic_word_terms!(Xnext, harmonic)
        for (word, coefficient) in residual_terms
          add_word_coefficient!(output_terms, word, weight * coefficient, counts)
        end
      end
      push!(wave, Xnext)
    end
  end

  return BlochWordPlan(support, order, wave, effective, counts)
end

function primitive_bloch_word_weight(word, inverse_weight, zero_harmonic)
  partial_sums = harmonic_partial_sums(word, zero_harmonic)
  isempty(partial_sums) && throw(ArgumentError("harmonic word must not be empty"))
  all(!iszero, partial_sums[1:(end - 1)]) ||
    throw(ArgumentError("primitive weight requires a first-return word"))
  iszero(last(partial_sums)) ||
    throw(ArgumentError("primitive weight requires a closed word"))

  weight = one(inverse_weight(first(partial_sums)))
  for partial_sum in partial_sums[1:(end - 1)]
    weight *= inverse_weight(partial_sum)
  end
  return weight
end

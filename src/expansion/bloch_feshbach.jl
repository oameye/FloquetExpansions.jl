struct LyndonBracketNode
  left::Int
  right::Int
end

struct LyndonTerm{C}
  node::Int
  coefficient::C
end

struct LieSeriesPlan{H,C}
  letters::Vector{H}
  brackets::Vector{LyndonBracketNode}
  connected_log::Vector{Dict{H,Vector{LyndonTerm{C}}}}
  effective::Vector{Vector{LyndonTerm{C}}}
end

struct LyndonBracketTable{H,C}
  letter_ids::Dict{H,Int}
  word_ids::Dict{Vector{H},Int}
  brackets::Vector{LyndonBracketNode}
  expansions::Dict{Vector{H},WordPolynomial{H,C}}
end

function LyndonBracketTable(letters::Vector{H}, ::Type{WordPolynomial{H,C}}) where {H,C}
  letter_ids = Dict{H,Int}(letter => index for (index, letter) in enumerate(letters))
  return LyndonBracketTable(
    letter_ids,
    Dict{Vector{H},Int}(),
    LyndonBracketNode[],
    Dict{Vector{H},WordPolynomial{H,C}}(),
  )
end

function lyndon_node!(table::LyndonBracketTable{H}, word::Vector{H}) where {H}
  length(word) == 1 && return table.letter_ids[only(word)]
  haskey(table.word_ids, word) && return table.word_ids[word]
  prefix, suffix = lyndon_standard_factorization(word)
  left = lyndon_node!(table, prefix)
  right = lyndon_node!(table, suffix)
  push!(table.brackets, LyndonBracketNode(left, right))
  node = length(table.letter_ids) + length(table.brackets)
  table.word_ids[word] = node
  return node
end

function lyndon_terms!(
  table::LyndonBracketTable{H,C}, polynomial::WordPolynomial{H,C}
) where {H,C}
  return LyndonTerm{C}[
    LyndonTerm(lyndon_node!(table, word), coefficient) for
    (word, coefficient) in lyndon_coordinates(polynomial, table.expansions)
  ]
end

function lyndon_series_terms!(
  table::LyndonBracketTable{H,C}, series::Dict{H,WordPolynomial{H,C}}
) where {H,C}
  compiled = Dict{H,Vector{LyndonTerm{C}}}()
  for (harmonic, polynomial) in series
    iszero(polynomial) || (compiled[harmonic] = lyndon_terms!(table, polynomial))
  end
  return compiled
end

function compile_lie_series_plan(support::Vector{Int}, order::Int)
  letters = sort!(unique(support))
  words = van_vleck_words(letters, order)
  table = LyndonBracketTable(letters, eltype(words.effective))
  connected_log = [lyndon_series_terms!(table, series) for series in words.connected_log]
  effective = [lyndon_terms!(table, polynomial) for polynomial in words.effective]
  return LieSeriesPlan(letters, table.brackets, connected_log, effective)
end

function evaluate_lyndon_brackets(
  plan::LieSeriesPlan{H}, components::AbstractDict{H,T}, product
) where {H,T}
  nletters = length(plan.letters)
  bracket_values = Vector{T}(undef, nletters + length(plan.brackets))
  for (index, letter) in enumerate(plan.letters)
    bracket_values[index] = components[letter]
  end
  for (index, bracket) in enumerate(plan.brackets)
    left = bracket_values[bracket.left]
    right = bracket_values[bracket.right]
    commutator = product(left, right) - product(right, left)
    bracket_values[nletters + index] = SQA.simplify(commutator)::T
  end
  return bracket_values
end

function lie_combination(
  terms::Vector{LyndonTerm{C}}, bracket_values::Vector{T}, zero_component::T, scale::Number
) where {C,T}
  value = zero_component
  for term in terms
    value += (scale * term.coefficient) * bracket_values[term.node]
  end
  return SQA.simplify(value)::T
end

function lie_series_component(
  output::Dict{H,Vector{LyndonTerm{C}}},
  bracket_values::Vector{T},
  zero_component::T,
  scale::Number,
) where {H,C,T}
  series = Dict{H,T}()
  for (harmonic, terms) in output
    component = lie_combination(terms, bracket_values, zero_component, scale)
    iszero(component) || (series[harmonic] = component)
  end
  return series
end

function evaluate_lie_series(
  plan::LieSeriesPlan{H},
  components::AbstractDict{H,T},
  zero_component::T,
  convention::ComponentConvention,
) where {H,T}
  (; product, generator_phase) = convention
  order_phase = im * generator_phase
  micromotion_phase = conj(generator_phase)
  bracket_values = evaluate_lyndon_brackets(plan, components, product)
  effective = T[
    lie_combination(terms, bracket_values, zero_component, order_phase^(n - 1)) for
    (n, terms) in enumerate(plan.effective)
  ]
  micromotion_series = [
    lie_series_component(
      output, bracket_values, zero_component, micromotion_phase * order_phase^n
    ) for (n, output) in enumerate(plan.connected_log)
  ]
  return (; micromotion_series, effective)
end

function van_vleck_expansion(
  generator::P, gauge::VanVleck{BlochFeshbach}, order::Int, provenance::R
) where {P<:PeriodicGenerator,R<:FloquetProvenance}
  (; components, zero_component, wd) = generator
  plan = compile_lie_series_plan(collect(keys(generator)), order)
  (; micromotion_series, effective) = evaluate_lie_series(
    plan, components, zero_component, component_convention(generator)
  )
  micromotion_components = P[
    periodic_generator(series, wd, zero_component)::P for series in micromotion_series
  ]

  return FloquetExpansion(
    generator, micromotion_components, effective, gauge, order, Uncompleted(), provenance
  )
end

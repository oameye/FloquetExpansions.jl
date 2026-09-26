struct BlochConventions{F,W<:Number,K<:Number}
  product::F
  weight_phase::W
  kick_phase::K
end

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

function compile_lie_series_plan(support::Vector{H}, order::Int, zero_harmonic::H) where {H}
  letters = sort!(unique(support))
  words = van_vleck_words(letters, order, zero_harmonic)
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
  conventions::BlochConventions,
) where {H,T}
  (; product, weight_phase, kick_phase) = conventions
  bracket_values = evaluate_lyndon_brackets(plan, components, product)
  effective = T[
    lie_combination(terms, bracket_values, zero_component, weight_phase^(n - 1)) for
    (n, terms) in enumerate(plan.effective)
  ]
  kick_series = [
    lie_series_component(
      output, bracket_values, zero_component, kick_phase * weight_phase^n
    ) for (n, output) in enumerate(plan.connected_log)
  ]
  return (; kick_series, effective)
end

bloch_conventions(::PeriodicGenerator{SQA.QAdd}) = BlochConventions(*, 1, im)
bloch_conventions(::PeriodicGenerator{Liouvillian}) = BlochConventions(compose, im, 1)

function floquet_expansion_impl(
  generator::PeriodicGenerator{SQA.QAdd},
  gauge::VanVleck{BlochFeshbach},
  order::Int,
  provenance::R,
) where {R<:FloquetProvenance}
  return bloch_feshbach_expansion(
    generator, gauge, order, provenance, bloch_conventions(generator)
  )
end

function floquet_expansion_impl(
  generator::PeriodicGenerator{Liouvillian},
  gauge::VanVleck{BlochFeshbach},
  order::Int,
  provenance::R,
) where {R<:FloquetProvenance}
  return bloch_feshbach_expansion(
    generator, gauge, order, provenance, bloch_conventions(generator)
  )
end

function bloch_feshbach_expansion(
  generator::P, gauge::G, order::Int, provenance::R, conventions::BlochConventions
) where {T,P<:PeriodicGenerator{T},G<:VanVleck{BlochFeshbach},R<:FloquetProvenance}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  generator isa PeriodicGenerator{SQA.QAdd} && require_hermitian_drive(generator)

  (; components, zero_component, wd) = generator
  plan = compile_lie_series_plan(collect(keys(generator)), order, 0)
  (; kick_series, effective) = evaluate_lie_series(
    plan, components, zero_component, conventions
  )
  kick_components = P[
    periodic_generator(series, wd, zero_component)::P for series in kick_series
  ]

  return FloquetExpansion(
    generator, kick_components, effective, gauge, order, Uncompleted(), provenance
  )
end

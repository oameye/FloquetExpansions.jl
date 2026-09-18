using FloquetExpansions

include(joinpath(@__DIR__, "connected_word_log.jl"))
include(joinpath(@__DIR__, "sparse_connected_word_log.jl"))

struct LyndonBracketNode
  left::Int
  right::Int
end

struct LyndonOutputTerm{C}
  node::Int
  coefficient::C
end

struct LyndonLogEvaluationPlan{H,C}
  leaves::Vector{H}
  brackets::Vector{LyndonBracketNode}
  outputs::Vector{Dict{H,Vector{LyndonOutputTerm{C}}}}
end

function lyndon_ensure_node!(
  word::Tuple,
  leaves::Vector{H},
  leaf_ids::Dict{H,Int},
  word_ids::Dict{Tuple,Int},
  brackets::Vector{LyndonBracketNode},
) where {H}
  if length(word) == 1
    return leaf_ids[first(word)]
  end
  haskey(word_ids, word) && return word_ids[word]
  left_word, right_word = lyndon_standard_factorization(word)
  left = lyndon_ensure_node!(left_word, leaves, leaf_ids, word_ids, brackets)
  right = lyndon_ensure_node!(right_word, leaves, leaf_ids, word_ids, brackets)
  push!(brackets, LyndonBracketNode(left, right))
  node = length(leaves) + length(brackets)
  word_ids[word] = node
  return node
end

function compile_lyndon_log_evaluation_plan(support, order::Int)
  C = Rational{Int}
  leaves = sort(unique(collect(support)))
  leaf_ids = Dict(harmonic => index for (index, harmonic) in enumerate(leaves))
  word_ids = Dict{Tuple,Int}()
  brackets = LyndonBracketNode[]
  bracket_cache = Dict{Tuple,HarmonicWordPolynomial{C}}()
  log_embedding = compile_sparse_connected_log_words(leaves, order)
  H = eltype(leaves)
  outputs = Vector{Dict{H,Vector{LyndonOutputTerm{C}}}}()

  for embedding in log_embedding
    compiled = Dict{H,Vector{LyndonOutputTerm{C}}}()
    for (harmonic, word_terms) in embedding
      polynomial = HarmonicWordPolynomial{C}(word_terms)
      terms = LyndonOutputTerm{C}[]
      for (word, coefficient) in lyndon_decomposition(polynomial, bracket_cache)
        node = lyndon_ensure_node!(word, leaves, leaf_ids, word_ids, brackets)
        push!(terms, LyndonOutputTerm(node, coefficient))
      end
      compiled[harmonic] = terms
    end
    push!(outputs, compiled)
  end

  return LyndonLogEvaluationPlan(leaves, brackets, outputs)
end

function evaluate_lyndon_log_plan(
  plan::LyndonLogEvaluationPlan{H,C},
  components::AbstractDict{H,T};
  product,
  zero_component::T,
  simplifier=identity,
) where {H,C,T}
  nleaves = length(plan.leaves)
  values = Vector{T}(undef, nleaves + length(plan.brackets))
  for (index, harmonic) in enumerate(plan.leaves)
    values[index] = components[harmonic]
  end

  for (index, bracket) in enumerate(plan.brackets)
    left = values[bracket.left]
    right = values[bracket.right]
    values[nleaves + index] = simplifier(product(left, right) - product(right, left))::T
  end

  output = Vector{Dict{H,T}}()
  for embedding in plan.outputs
    evaluated = Dict{H,T}()
    for (harmonic, terms) in embedding
      value = zero_component
      for term in terms
        value += term.coefficient * values[term.node]
      end
      component = simplifier(value)::T
      iszero(component) || (evaluated[harmonic] = component)
    end
    push!(output, evaluated)
  end
  return output
end

lyndon_backend_products(plan::LyndonLogEvaluationPlan) = 2 * length(plan.brackets)

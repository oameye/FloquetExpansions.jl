struct CompiledHarmonicWord{H}
  harmonics::Tuple{Vararg{H}}
  partial_sums::Tuple{Vararg{H}}
  return_positions::Tuple{Vararg{Int}}
end

struct HarmonicWordPlanCounts
  raw_words::Int
  attempted_edges::Int
  pruned_edges::Int
  closed_words::Int
  first_return_words::Int
  folded_words::Int
  dag_nodes::Int
  dag_edges::Int
end

struct HarmonicWordPlan{H}
  support::Vector{H}
  word_length::Int
  closed::Vector{CompiledHarmonicWord{H}}
  first_return_indices::Vector{Int}
  folded_indices::Vector{Int}
  dag_nodes::Set{Tuple{Int,H}}
  dag_edges::Set{Tuple{Int,H,H}}
  counts::HarmonicWordPlanCounts
end

function harmonic_partial_sums(harmonics, zero_harmonic)
  sums = typeof(zero_harmonic)[]
  total = zero_harmonic
  for harmonic in harmonics
    total += harmonic
    push!(sums, total)
  end
  return sums
end

function compiled_harmonic_word(harmonics, zero_harmonic)
  sums = harmonic_partial_sums(harmonics, zero_harmonic)
  returns = findall(iszero, sums)
  H = typeof(zero_harmonic)
  return CompiledHarmonicWord{H}(Tuple(harmonics), Tuple(sums), Tuple(returns))
end

function is_closed(word::CompiledHarmonicWord)
  return !isempty(word.partial_sums) && iszero(last(word.partial_sums))
end
function is_first_return(word::CompiledHarmonicWord)
  is_closed(word) || return false
  return all(!iszero(sum) for sum in word.partial_sums[1:(end - 1)])
end

function fold_positions(word::CompiledHarmonicWord)
  is_closed(word) || return ()
  return Tuple(
    position for position in word.return_positions if position < length(word.harmonics)
  )
end

function return_blocks(word::CompiledHarmonicWord)
  is_closed(word) || throw(ArgumentError("return blocks require a closed harmonic word"))
  blocks = Tuple{Vararg{eltype(word.harmonics)}}[]
  first_index = 1
  for last_index in word.return_positions
    push!(blocks, word.harmonics[first_index:last_index])
    first_index = last_index + 1
  end
  return blocks
end

function reachable_harmonic_sums(
  support::Vector{H}, max_steps::Int, zero_harmonic::H
) where {H}
  reachable = [Set{H}() for _ in 0:max_steps]
  push!(reachable[1], zero_harmonic)
  for steps in 1:max_steps
    for partial_sum in reachable[steps], harmonic in support
      push!(reachable[steps + 1], partial_sum + harmonic)
    end
  end
  return reachable
end

function compile_harmonic_word_plan(
  support_input, word_length::Int; zero_harmonic=zero(first(support_input))
)
  word_length >= 1 || throw(ArgumentError("word length must be >= 1"))
  isempty(support_input) && throw(ArgumentError("harmonic support must not be empty"))

  support = unique(collect(support_input))
  H = typeof(zero_harmonic)
  all(harmonic isa H for harmonic in support) ||
    throw(ArgumentError("harmonic support and zero harmonic must have one common type"))
  support = H[harmonic for harmonic in support]

  reachable = reachable_harmonic_sums(support, word_length, zero_harmonic)
  words = CompiledHarmonicWord{H}[]
  first_return_indices = Int[]
  folded_indices = Int[]
  dag_nodes = Set{Tuple{Int,H}}([(0, zero_harmonic)])
  dag_edges = Set{Tuple{Int,H,H}}()
  harmonics = H[]
  partial_sums = H[]
  attempted_edges = Ref(0)
  pruned_edges = Ref(0)

  function visit(depth::Int, current_sum)
    for harmonic in support
      attempted_edges[] += 1
      next_sum = current_sum + harmonic
      next_depth = depth + 1
      remaining = word_length - next_depth
      can_close = if iszero(remaining)
        iszero(next_sum)
      else
        -next_sum in reachable[remaining + 1]
      end
      if !can_close
        pruned_edges[] += 1
        continue
      end

      push!(dag_nodes, (next_depth, next_sum))
      push!(dag_edges, (depth, current_sum, harmonic))
      push!(harmonics, harmonic)
      push!(partial_sums, next_sum)

      if next_depth == word_length
        returns = findall(iszero, partial_sums)
        word = CompiledHarmonicWord{H}(
          Tuple(harmonics), Tuple(partial_sums), Tuple(returns)
        )
        push!(words, word)
        if is_first_return(word)
          push!(first_return_indices, length(words))
        else
          push!(folded_indices, length(words))
        end
      else
        visit(next_depth, next_sum)
      end

      pop!(partial_sums)
      pop!(harmonics)
    end
    return nothing
  end

  visit(0, zero_harmonic)
  counts = HarmonicWordPlanCounts(
    length(support)^word_length,
    attempted_edges[],
    pruned_edges[],
    length(words),
    length(first_return_indices),
    length(folded_indices),
    length(dag_nodes),
    length(dag_edges),
  )
  return HarmonicWordPlan(
    support,
    word_length,
    words,
    first_return_indices,
    folded_indices,
    dag_nodes,
    dag_edges,
    counts,
  )
end

first_return_words(plan::HarmonicWordPlan) = plan.closed[plan.first_return_indices]
folded_words(plan::HarmonicWordPlan) = plan.closed[plan.folded_indices]

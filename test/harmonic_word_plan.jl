using Test

include(joinpath(@__DIR__, "helpers", "harmonic_word_plan.jl"))

function brute_closed_words(support, word_length, zero_harmonic)
  factors = ntuple(_ -> support, word_length)
  closed = Set{Tuple}()
  for harmonics in Iterators.product(factors...)
    partial_sums = harmonic_partial_sums(harmonics, zero_harmonic)
    iszero(last(partial_sums)) && push!(closed, Tuple(harmonics))
  end
  return closed
end

function planned_words(words)
  return Set{Tuple}(word.harmonics for word in words)
end

@testset "harmonic word plan closes and partitions exactly" begin
  support = [-1, 0, 1]
  for word_length in 1:6
    plan = compile_harmonic_word_plan(support, word_length)
    brute = brute_closed_words(support, word_length, 0)

    @test planned_words(plan.closed) == brute
    @test plan.counts.raw_words == length(support)^word_length
    @test plan.counts.closed_words == length(brute)
    @test plan.counts.first_return_words + plan.counts.folded_words == length(brute)
    @test isempty(intersect(planned_words(first_return_words(plan)), planned_words(folded_words(plan))))
    @test union(planned_words(first_return_words(plan)), planned_words(folded_words(plan))) == brute
    @test all(is_first_return, first_return_words(plan))
    @test all(!is_first_return(word) for word in folded_words(plan))
  end
end

@testset "folds are unique model-space return decompositions" begin
  plan = compile_harmonic_word_plan([-2, -1, 0, 1, 2], 5)
  for word in plan.closed
    blocks = return_blocks(word)
    @test Tuple(Iterators.flatten(blocks)) == word.harmonics
    for block in blocks
      block_word = compiled_harmonic_word(block, 0)
      @test is_closed(block_word)
      if length(block) == 1
        @test iszero(only(block))
      else
        @test is_first_return(block_word)
      end
    end
    @test isempty(fold_positions(word)) == is_first_return(word)
  end
end

@testset "first-return words reproduce Bloch bare-excursion support" begin
  support = [-2, -1, 0, 1, 2]
  nonzero_support = filter(!iszero, support)

  order_two = compile_harmonic_word_plan(support, 2)
  expected_B1 = Set((harmonic, -harmonic) for harmonic in nonzero_support)
  @test planned_words(first_return_words(order_two)) == expected_B1

  order_three = compile_harmonic_word_plan(support, 3)
  expected_B2_primitive = Set{Tuple}()
  expected_B2_fold = Set{Tuple}()
  for outer_harmonic in nonzero_support, inner_harmonic in nonzero_support
    middle_harmonic = outer_harmonic - inner_harmonic
    middle_harmonic in support || continue
    push!(
      expected_B2_primitive,
      (inner_harmonic, middle_harmonic, -outer_harmonic),
    )
  end
  for harmonic in nonzero_support
    push!(expected_B2_fold, (0, harmonic, -harmonic))
  end

  # Flow words are written in propagation order, i.e. opposite to the operator-product order.
  # Thus the primitive word (k, m-k, -m) represents H[-m] H[m-k] H[k], exactly the
  # first term of the certified B₂ Bloch coefficient in #96.
  @test planned_words(first_return_words(order_three)) == expected_B2_primitive
  @test issubset(expected_B2_fold, planned_words(folded_words(order_three)))
end

@testset "support pruning and frequency-state DAG are exact" begin
  support = [-2, 0, 1]
  plan = compile_harmonic_word_plan(support, 6)
  brute = brute_closed_words(support, 6, 0)
  full_prefix_tree_nodes = sum(length(support)^depth for depth in 0:6)
  full_prefix_tree_edges = full_prefix_tree_nodes - 1

  @test planned_words(plan.closed) == brute
  @test plan.counts.pruned_edges > 0
  @test plan.counts.dag_nodes < full_prefix_tree_nodes
  @test plan.counts.dag_edges < full_prefix_tree_edges
  @test length(plan.dag_nodes) == plan.counts.dag_nodes
  @test length(plan.dag_edges) == plan.counts.dag_edges

  permuted = compile_harmonic_word_plan(reverse(support), 6)
  @test planned_words(permuted.closed) == planned_words(plan.closed)
  @test planned_words(first_return_words(permuted)) == planned_words(first_return_words(plan))
  @test permuted.counts.dag_nodes == plan.counts.dag_nodes
  @test permuted.counts.dag_edges == plan.counts.dag_edges
end

struct TestHarmonicIndex
  n1::Int
  n2::Int
end

Base.zero(::Type{TestHarmonicIndex}) = TestHarmonicIndex(0, 0)
Base.zero(::TestHarmonicIndex) = zero(TestHarmonicIndex)
Base.iszero(index::TestHarmonicIndex) = iszero(index.n1) && iszero(index.n2)
function Base.:+(left::TestHarmonicIndex, right::TestHarmonicIndex)
  return TestHarmonicIndex(left.n1 + right.n1, left.n2 + right.n2)
end
Base.:-(index::TestHarmonicIndex) = TestHarmonicIndex(-index.n1, -index.n2)

@testset "harmonic word plan is not scalar-frequency specific" begin
  zero_index = zero(TestHarmonicIndex)
  e1 = TestHarmonicIndex(1, 0)
  e2 = TestHarmonicIndex(0, 1)
  support = [zero_index, e1, -e1, e2, -e2]
  plan = compile_harmonic_word_plan(support, 4; zero_harmonic=zero_index)
  brute = brute_closed_words(support, 4, zero_index)

  @test planned_words(plan.closed) == brute
  @test all(iszero(last(word.partial_sums)) for word in plan.closed)
  @test plan.counts.first_return_words > 0
  @test plan.counts.folded_words > 0
end

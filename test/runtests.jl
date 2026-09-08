using FloquetExpansions
using ParallelTestRunner: ParallelTestRunner

testsuite = ParallelTestRunner.find_tests(@__DIR__)
args = ParallelTestRunner.parse_args(ARGS)

if ParallelTestRunner.filter_tests!(testsuite, args)
  delete!(testsuite, "quality/JET")
end

# Expert API is stable but intentionally not exported. Import it explicitly into each isolated
# test sandbox so tests exercise the qualified public surface without widening package exports.
init_code = quote
  using FloquetExpansions: Completion,
    Uncompleted,
    CompletionAlgorithm,
    CompletionFactorization,
    GramStage,
    GramFactorization,
    SpectralFactorization,
    CompletionObstruction,
    FractionalJumpOnset
end

ParallelTestRunner.runtests(FloquetExpansions, args; testsuite, init_code)

using Test
using FloquetExpansions

include(joinpath(@__DIR__, "helpers", "connected_word_log.jl"))
include(joinpath(@__DIR__, "helpers", "sparse_connected_word_log.jl"))

@testset "sparse connected log matches generic word oracle" begin
  support = [-1, 0, 1]
  order = 6
  _, _, generic = connected_word_reconstruction(support, order)
  sparse = compile_sparse_connected_log_words(support, order)

  @test length(sparse) == length(generic.log_embedding)
  for n in eachindex(sparse)
    harmonics = union(keys(sparse[n]), keys(generic.log_embedding[n]))
    for harmonic in harmonics
      sparse_terms = get(sparse[n], harmonic, Dict{Tuple,Rational{Int}}())
      generic_terms = get(
        generic.log_embedding[n],
        harmonic,
        HarmonicWordPolynomial{Rational{Int}}(),
      ).terms
      @test sparse_terms == generic_terms
    end
  end
end

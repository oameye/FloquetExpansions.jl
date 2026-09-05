module FloquetExpansions

using Reexport: @reexport
using LinearAlgebra: LinearAlgebra
using Symbolics: Symbolics

using SecondQuantizedAlgebra: SecondQuantizedAlgebra
@reexport using SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra

# `@public` in SQA but NOT exported, so `@reexport` does not forward them.
using SecondQuantizedAlgebra: expim, exponential_form, trigonometric_form
export expim, exponential_form, trigonometric_form

include("periodic_operator.jl")
include("completion_types.jl")
include("matrix_series.jl")
include("liouvillian.jl")
include("quasienergy.jl")
include("engine.jl")
include("gksl_coordinates.jl")
include("gram_completion.jl")
include("gram_recursion.jl")
include("completion.jl")
include("gksl_floquet.jl")

export PeriodicGenerator, Gauge, VanVleck, QuasienergyOperator, harmonic_range
export time_average, derivative, antiderivative, support, harmonics
export FloquetExpansion,
  floquet_expansion, effective_generator, effective_component, micromotion
export Completion, Uncompleted, CompletionAlgorithm, Gram, Spectral, CompletionFactorization
export GramFactorization, CompletionObstruction, FractionalJumpOnset
export positive_completion,
  dissipative_frame, channels, positivity_conditions, regularity_conditions, factorization
export Liouvillian,
  liouvillian, terms, hamiltonian_action, dissipator, compose, collapse, jump
export DissipativeFrame,
  hamiltonian, hamiltonian_component, kossakowski, kossakowski_component

end # module FloquetExpansions

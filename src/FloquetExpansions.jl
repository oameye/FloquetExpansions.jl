module FloquetExpansions

using Reexport: @reexport
using LinearAlgebra: LinearAlgebra
using SciMLPublic: @public
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
include("completion_linear_algebra.jl")
include("liouvillian.jl")
include("quasienergy.jl")
include("engine.jl")
include("gksl_coordinates.jl")
include("completion_conversion.jl")
include("completion_frame.jl")
include("gram_completion.jl")
include("gram_recursion.jl")
include("spectral_completion.jl")
include("completion.jl")
include("gksl_floquet.jl")

export PeriodicGenerator, Gauge, VanVleck, QuasienergyOperator, harmonic_range
export time_average, derivative, antiderivative, support, harmonics
export FloquetExpansion,
  floquet_expansion, order, effective_generator, effective_component, micromotion
export Gram, Spectral
export positive_completion,
  dissipative_frame, channels, positivity_conditions, regularity_conditions, factorization
export Liouvillian,
  liouvillian, terms, hamiltonian_action, dissipator, compose, collapse, jump
export DissipativeFrame,
  hamiltonian, hamiltonian_component, kossakowski, kossakowski_component

# Stable expert API that is intentionally qualified rather than exported.
@public Completion,
  Uncompleted,
  CompletionAlgorithm,
  CompletionFactorization,
  GramStage,
  GramFactorization,
  SpectralFactorization,
  CompletionObstruction,
  FractionalJumpOnset

end # module FloquetExpansions

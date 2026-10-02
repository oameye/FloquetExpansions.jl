module FloquetExpansions

using Reexport: @reexport
using LinearAlgebra: LinearAlgebra
using SciMLPublic: @public
using Symbolics: Symbolics

using SecondQuantizedAlgebra: SecondQuantizedAlgebra
@reexport using SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra

using SecondQuantizedAlgebra: expim, exponential_form, trigonometric_form
export expim, exponential_form, trigonometric_form

include("generators/periodic_generator.jl")
include("generators/liouvillian.jl")
include("generators/channels.jl")
include("generators/quasienergy.jl")

include("words/lyndon_words.jl")
include("words/van_vleck_words.jl")

include("expansion/gauges.jl")
include("expansion/floquet_expansion.jl")
include("expansion/hori_deprit.jl")
include("expansion/bloch_feshbach.jl")

include("gksl/dissipative_frame.jl")
include("gksl/coordinates.jl")
include("gksl/floquet.jl")
include("gksl/gram_kernel.jl")

include("native/graded_channels.jl")
include("native/static_solve.jl")
include("native/exact_static_solve.jl")
include("native/static_step.jl")
include("native/dense_representation.jl")
include("native/exact_representation.jl")
include("native/charge_grading.jl")
include("native/operator_algebra.jl")
include("native/algebraic_representation.jl")
include("native/sqa_lowering.jl")
include("native/recurrence.jl")
include("native/driver.jl")
include("native/symbolic_parameters.jl")
include("native/realization.jl")

include("completion/types.jl")
include("completion/backend/matrix_series.jl")
include("completion/backend/linear_algebra.jl")
include("completion/backend/conversion.jl")
include("completion/frame_discovery.jl")
include("completion/accessors.jl")
include("completion/gram/factorization.jl")
include("completion/gram/recursion.jl")
include("completion/spectral.jl")
include("completion/positive_completion.jl")

export BlochFeshbach, HoriDeprit
export GKSLNormalForm,
  Gauge, PeriodicGenerator, QuasienergyOperator, VanVleck, harmonic_range
export antiderivative, derivative, harmonics, support, time_average
export FloquetExpansion,
  effective_component, effective_generator, floquet_expansion, micromotion, order
export Gram, Spectral
export channels,
  dissipative_frame,
  factorization,
  positive_completion,
  positivity_conditions,
  regularity_conditions
export Liouvillian,
  collapse, compose, dissipator, hamiltonian_action, jump, liouvillian, terms
export DissipativeFrame,
  hamiltonian, hamiltonian_component, kossakowski, kossakowski_component

@public Completion,
CompletionAlgorithm,
CompletionFactorization,
GramFactorization,
GramStage,
SpectralFactorization,
Uncompleted,
ExpansionAlgorithm

end

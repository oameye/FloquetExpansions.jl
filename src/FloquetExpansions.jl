module FloquetExpansions

using Reexport: @reexport
using LinearAlgebra: LinearAlgebra
using Nemo: Nemo
using SciMLPublic: @public
using Symbolics: Symbolics

using SecondQuantizedAlgebra: SecondQuantizedAlgebra
@reexport using SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra

using SecondQuantizedAlgebra: expim, exponential_form, trigonometric_form
export expim, exponential_form, trigonometric_form
using SecondQuantizedAlgebra: channels, hamiltonian

include("generators/periodic_generator.jl")
include("generators/superoperator.jl")
include("generators/channels.jl")
include("generators/quasienergy.jl")

include("words/lyndon_words.jl")
include("words/van_vleck_words.jl")

include("expansion/gauges.jl")
include("expansion/floquet_expansion.jl")
include("expansion/hori_deprit.jl")
include("expansion/bloch_feshbach.jl")

include("gksl/coefficients.jl")
include("gksl/floquet.jl")

include("native/graded_channels.jl")
include("native/exact_static_solve.jl")
include("native/static_step.jl")
include("native/operator_algebra.jl")
include("native/algebraic_representation.jl")
include("native/sqa_lowering.jl")
include("native/recurrence.jl")
include("native/driver.jl")
include("native/parameter_field.jl")
include("native/realization.jl")

include("kraus/time_kernel.jl")
include("kraus/no_jump_frame.jl")
include("kraus/period_map.jl")

include("completion/types.jl")
include("completion/backend/matrix_series.jl")
include("completion/backend/linear_algebra.jl")
include("completion/backend/conversion.jl")
include("completion/basis_discovery.jl")
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
export PeriodMap, PeriodMapExpansion, period_map, superoperator, trace_defect
export JumpAmplitude, amplitude, jump_amplitudes, jump_times, weight
export Gram, Spectral
export channels,
  jump_basis,
  factorization,
  positive_completion,
  positivity_conditions,
  regularity_conditions
export hamiltonian, hamiltonian_component, kossakowski_component

@public Completion,
CompletionAlgorithm,
CompletionFactorization,
CompletionObstruction,
FractionalJumpOnset,
GramFactorization,
GramStage,
SpectralFactorization,
Uncompleted,
ExpansionAlgorithm

end

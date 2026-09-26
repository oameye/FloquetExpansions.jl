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

include("gauges.jl")
include("periodic_operator.jl")
include("completion_types.jl")
include("matrix_series.jl")
include("completion_linear_algebra.jl")
include("liouvillian.jl")
include("quasienergy.jl")
include("floquet_expansion.jl")
include("hori_deprit.jl")
include("bloch_feshbach/lyndon_words.jl")
include("bloch_feshbach/van_vleck_words.jl")
include("bloch_feshbach/bloch_feshbach.jl")
include("gksl_coordinates.jl")
include("completion_conversion.jl")
include("completion_frame.jl")
include("gram_completion.jl")
include("gram_recursion.jl")
include("spectral_completion.jl")
include("completion.jl")
include("gksl_floquet.jl")

export BlochFeshbach, HoriDeprit
export Gauge, PeriodicGenerator, QuasienergyOperator, VanVleck, harmonic_range
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
Uncompleted, ExpansionAlgorithm

end

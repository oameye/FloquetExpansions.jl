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
include("liouvillian.jl")
include("quasienergy.jl")
include("engine.jl")

export PeriodicGenerator, Gauge, VanVleck, QuasienergyOperator, harmonic_range
export time_average, derivative, antiderivative, support, harmonics
export FloquetExpansion, floquet_expansion, effective_generator, micromotion
export effective_hamiltonian, kick_operator
export Liouvillian, actions, hamiltonian_action, dissipator, compose, collapse, jump

end # module FloquetExpansions

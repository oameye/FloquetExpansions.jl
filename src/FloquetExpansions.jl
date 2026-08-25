module FloquetExpansions

using Reexport: @reexport
using LinearAlgebra: LinearAlgebra
using Symbolics: Symbolics

# Second quantized algebra
import SecondQuantizedAlgebra as SQA
@reexport using SecondQuantizedAlgebra

# `@public` in SQA but NOT exported, so `@reexport` does not forward them. Users need `expim` and
# `exponential_form` to write a drive, and to read `kick_operator(vv, t)` back.
using SecondQuantizedAlgebra: expim, exponential_form, trigonometric_form
export expim, exponential_form, trigonometric_form

include("periodic_operator.jl")
include("quasienergy.jl")
include("collector.jl")
include("engine.jl")

export PeriodicOperator, Gauge, VanVleck, QuasienergyOperator, harmonic_range
export time_average, derivative, antiderivative, support, harmonics
export FloquetExpansion, floquet_expansion, effective_hamiltonian, kick_operator

end # module FloquetExpansions

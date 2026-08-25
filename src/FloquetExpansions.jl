module FloquetExpansions

using Reexport: @reexport
using LinearAlgebra: LinearAlgebra
using Symbolics: Symbolics

# Second quantized algebra
import SecondQuantizedAlgebra as SQA
@reexport using SecondQuantizedAlgebra

include("periodic_operator.jl")
include("collector.jl")
include("engine.jl")

export PeriodicOperator, Gauge, VanVleck
export time_average, derivative, antiderivative, support, harmonics
export FloquetExpansion, floquet_expansion, effective_hamiltonian, kick_operator

end # module FloquetExpansions

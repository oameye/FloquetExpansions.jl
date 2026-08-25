module FloquetExpansions

using Reexport: @reexport
using LinearAlgebra: LinearAlgebra

# Second quantized algebra
import SecondQuantizedAlgebra as SQA
@reexport using SecondQuantizedAlgebra

include("periodic_operator.jl")

export PeriodicOperator, Gauge, VanVleck
export time_average, derivative, antiderivative, support

end # module FloquetExpansions

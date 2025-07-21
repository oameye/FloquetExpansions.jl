"""
$(DocStringExtensions.README)
"""
module FloquetExpansions
using DocStringExtensions: DocStringExtensions
using Reexport: @reexport

# symbolic utilities
using Symbolics: Symbolics
using SymbolicUtils: SymbolicUtils

# trigometric symbolic utilities
using QuestBase: QuestBase

# Second quantized algebra
using SecondQuantizedAlgebra:
    SecondQuantizedAlgebra, QMul, QAdd, QTerm, Create, Destroy, QSym
import SecondQuantizedAlgebra as SQA

# Van Vleck recursion formula
using VanVleckRecursion: VanVleckRecursion

include("utils.jl")
include("rotate.jl")
include("fourier_components.jl")
include("frequency_expansion.jl")

@reexport using SecondQuantizedAlgebra
export rotating_wave_approximation
export quesienergy_operator
export rotate

end # module FloquetExpansions

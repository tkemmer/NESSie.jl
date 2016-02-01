module Nonlocal

using ..ProteinES
import Distances: euclidean

export
    # bem.jl
    φmol, ∂ₙφmol, cauchy,
    # fem.jl
    espotential

# Include module files
include("nonlocal/bem.jl")
include("nonlocal/fem.jl")

end # module

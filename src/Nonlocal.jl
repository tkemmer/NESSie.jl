module Nonlocal

using ..ProteinES
import ..ProteinES: Radon, Rjasanow, eye!, gemv!, axpy!
import ..ProteinES.Local: φmol, ∂ₙφmol
import Distances: euclidean

export
    # bem.jl
    cauchy,
    # fem.jl
    espotential

# Include module files
include("nonlocal/bem.jl")
include("nonlocal/fem.jl")

end # module

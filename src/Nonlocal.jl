module Nonlocal

using ..ProteinES
using ..ProteinES: Radon, Rjasanow, eye!, gemv!, axpy!, reverseindex
using ..ProteinES.Local: φmol, ∂ₙφmol
using Distances: euclidean

export
    # bem.jl
    cauchy,
    # fem.jl
    espotential

# Include module files
include("nonlocal/bem.jl")
include("nonlocal/fem.jl")

end # module

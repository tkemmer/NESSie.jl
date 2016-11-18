module BEM

using ..ProteinES
using ..ProteinES: Radon, Rjasanow, eye!, pluseye!, gemv!, axpy!, reverseindex, φmol, ∂ₙφmol
using Distances: euclidean

export
    # nonlocal.jl
    cauchy, φΩ

include("bem/local.jl")
include("bem/nonlocal.jl")

end # module

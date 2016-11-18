module BEM

using ..ProteinES
using ..ProteinES: Radon, Rjasanow, pluseye!, gemv!, gemv, gemm, axpy!, φmol, ∂ₙφmol
using Distances: euclidean

export solvelocal, solvenonlocal

abstract BEMResult{T <: AbstractFloat}

include("bem/local.jl")
include("bem/nonlocal.jl")

end # module

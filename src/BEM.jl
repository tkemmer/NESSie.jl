module BEM

using ..ProteinES
using ..ProteinES: Radon, Rjasanow, pluseye!, gemv!, gemv, gemm, axpy!, φmol, ∂ₙφmol
using Distances: euclidean

export solvelocal, solvenonlocal

abstract BEMResult{T <: AbstractFloat}

include("bem/local/solver.jl")
include("bem/local/potentials.jl")
include("bem/nonlocal/solver.jl")
include("bem/nonlocal/potentials.jl")

end # module

module BEM

using ..ProteinES
using ..ProteinES: Radon, Rjasanow, σ, potprefactor, pluseye!, gemv!, gemv, gemm, axpy!, yukawa
using Distances: euclidean

import ..ProteinES: φΩ, φΣ
export solve, φΩ, φΣ

abstract BEMResult{T <: AbstractFloat}

include("bem/local/solver.jl")
include("bem/local/potentials.jl")
include("bem/nonlocal/solver.jl")
include("bem/nonlocal/potentials.jl")

end # module

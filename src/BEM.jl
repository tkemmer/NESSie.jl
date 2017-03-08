module BEM

using ..ProteinES
using ..ProteinES: Radon, Rjasanow, σ, ec, potprefactor, pluseye!, gemv!, gemv, gemm, axpy!, yukawa
using Distances: euclidean

import ..ProteinES: φΩ, φΣ
export solve, φΩ, φΣ, rfenergy

abstract type BEMResult{T <: AbstractFloat} end

# local electrostatics
include("bem/local/solver.jl")
include("bem/local/potentials.jl")
include("bem/local/energy.jl")

# nonlocal electrostatics
include("bem/nonlocal/solver.jl")
include("bem/nonlocal/potentials.jl")

end # module

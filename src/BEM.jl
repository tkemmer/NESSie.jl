module BEM

using ..NESSie
using ..NESSie: σ, ec, potprefactor, pluseye!, gemv!, gemv, gemm, axpy!, yukawa
using Distances: euclidean
using LinearAlgebra: ⋅, rmul!

export solve, φΩ, φΣ, rfenergy

"""
    abstract type BEMResult{T, E <: SurfaceElement{T}} end

Abstract base type for all BEM solver results
"""
abstract type BEMResult{T, E <: SurfaceElement{T}} end

# local electrostatics
include("bem/local/solver.jl")
include("bem/local/potentials.jl")
include("bem/local/energy.jl")

# nonlocal electrostatics
include("bem/nonlocal/solver.jl")
include("bem/nonlocal/potentials.jl")
include("bem/nonlocal/energy.jl")

end # module

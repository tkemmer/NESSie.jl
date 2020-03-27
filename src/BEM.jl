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

include("bem/local.jl")
include("bem/nonlocal.jl")
include("bem/post.jl")

end # module

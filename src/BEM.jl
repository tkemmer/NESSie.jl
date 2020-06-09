module BEM

using ..NESSie
using ..NESSie: _etol, σ, ec, potprefactor, pluseye!, gemv!, gemv, gemm, axpy!, yukawa
using Distances: euclidean
using LinearAlgebra: ⋅, rmul!

export BEMResult, LocalBEMResult, NonlocalBEMResult, rfenergy, solve, φΩ, φΣ

"""
    abstract type BEMResult{T, E <: SurfaceElement{T}} end

Abstract base type for all BEM solver results
"""
abstract type BEMResult{T, E <: SurfaceElement{T}} end

include("bem/local.jl")
include("bem/nonlocal.jl")
include("bem/post.jl")

end # module

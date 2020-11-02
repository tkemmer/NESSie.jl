module BEM

using ..NESSie
using ..NESSie: _etol, σ, ec, potprefactor, pluseye!, gemv!, gemv, gemm, axpy!, yukawa
using Base: Threads
using Distances: euclidean
using ImplicitArrays: BlockMatrix, FixedValueArray, InteractionFunction, InteractionMatrix
using IterativeSolvers: gmres
using LinearAlgebra
using Preconditioners: DiagonalPreconditioner

export BEMResult, LocalBEMResult, NonlocalBEMResult, rfenergy, solve, solve_implicit, φΩ, φΣ

"""
    abstract type BEMResult{T, E <: SurfaceElement{T}} end

Abstract base type for all BEM solver results
"""
abstract type BEMResult{T, E <: SurfaceElement{T}} end

include("bem/implicit.jl")
include("bem/local.jl")
include("bem/nonlocal.jl")
include("bem/post.jl")

end # module

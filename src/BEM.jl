module BEM

using ..NESSie
using ..NESSie: _axpy!, _closest_element_id, _etol, _gemm, _gemv, _gemv!, _molpotential,
    _molpotential_dn, _pluseye!, σ, ec, potprefactor, yukawa
using AutoHashEquals
using Base: Threads
using Distances: euclidean
using ImplicitArrays: BlockMatrix, FixedValueArray, InteractionFunction, InteractionMatrix
using IterativeSolvers: gmres
using LinearAlgebra
using Preconditioners: DiagonalPreconditioner

export BEMResult, LocalBEMResult, NonlocalBEMResult, solve

"""
    abstract type BEMResult{T, E <: SurfaceElement{T}} end

Abstract base type for all BEM solver results
"""
abstract type BEMResult{T, E <: SurfaceElement{T}} end

@inline function Base.show(io::IO, ::MIME"text/plain", bem::BEMResult)
    show(io, bem)
end

@inline function Base.show(io::IO, bem::BEMResult)
    print(io,
        "$(typeof(bem))",
        "(nodes = ", length(bem.model.nodes),
        ", elements = ", length(bem.model.elements),
        ", charges = ", length(bem.model.charges), ")"
    )
end

include("bem/implicit.jl")
include("bem/local.jl")
include("bem/nonlocal.jl")
include("bem/post.jl")

include("bem/deprecation.jl")

end # module

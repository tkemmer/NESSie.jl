module FEM

using ..NESSie
using ..NESSie:  ec, φmol, ∇φmol, reverseindex, axpy!, potprefactor
using Distances: euclidean

export solve

"""
    abstract type FEMResult{T, E <: VolumeElement{T}} end

Abstract base type for all FEM solver results
"""
abstract type FEMResult{T, E <: VolumeElement{T}} end

# nonlocal electrostatics
include("fem/nonlocal/solver.jl")

end # module

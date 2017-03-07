module FEM

using ..NESSie
using ..NESSie:  ec, φmol, ∇φmol, reverseindex, axpy!, potprefactor
using Distances: euclidean

export espotential

include("fem/nonlocal.jl")

end # module

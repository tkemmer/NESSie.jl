module TestModel

using ..NESSie
using ..NESSie: _norm, ε0, ec, cos, potprefactor
using Distances: euclidean
using LinearAlgebra: norm, rmul!
using SpecialFunctions: besseli, besselk

#=
    Born models
=#
include("testmodel/born/model.jl")
export BornIon, bornion

include("testmodel/born/potentials.jl")
export φΩ, φΣ

#=
    Multi-charge spheres
=#
include("testmodel/xie/model.jl")
export XieModel

include("testmodel/xie/common.jl")

include("testmodel/xie/nonlocal1.jl")
export NonlocalXieModel1

include("testmodel/xie/potentials.jl")

end # module

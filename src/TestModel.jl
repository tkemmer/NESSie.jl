module TestModel

using ..NESSie
using ..NESSie: potprefactor, cos
using Distances: euclidean
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

end # module

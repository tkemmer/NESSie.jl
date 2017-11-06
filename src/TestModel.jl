module TestModel

using ..NESSie
using ..NESSie: potprefactor
using Distances: euclidean

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

end # module

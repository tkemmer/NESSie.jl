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

end # module

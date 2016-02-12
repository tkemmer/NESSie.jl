module Local

using ..ProteinES
using ..ProteinES: ddot
using Distances: euclidean

export φmol, ∂ₙφmol

include("local/misc.jl")
include("local/bem.jl")

end # module

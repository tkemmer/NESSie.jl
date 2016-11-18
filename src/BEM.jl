module BEM

using ..ProteinES
using ..ProteinES: Radon, Rjasanow, pluseye!, gemv!, gemv, gemm, axpy!, φmol, ∂ₙφmol
using Distances: euclidean

export
    # local.jl
    solvelocal,
    # nonlocal.jl
    cauchy, φΩ

include("bem/local.jl")
include("bem/nonlocal.jl")

end # module

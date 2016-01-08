module NonLocalBEM

import Base: cos, sign
import Base.LinAlg.BLAS: gemv!, axpy!
import JSON: json

export
    # types.jl
    Element, Triangle, Tetrahedron, Charge, Option, SingleLayer, DoubleLayer,
    # util.jl
    props!, indexmap, unpack, vertexnormals, xml3d_mesh, eye!, isdegenerate, cos, cathetus, sign, distance, reverseindex,
    # input file readers
    readhmo, readmatlab, readoff, readpqr,
    # bem.jl
    defaultopt, singularpot, cauchy

# Include module files
include("types.jl")
include("util.jl")
include("input/hmo.jl")
include("input/matlab.jl")
include("input/off.jl")
include("input/pqr.jl")
include("bem.jl")

# Default options
const defaultopt64 = Option(2., 78., 1.8, 20.)
const defaultopt32 = Option(2f0, 78f0, 1.8f0, 20f0)
defaultopt(::Type{Float64}) = defaultopt64
defaultopt(::Type{Float32}) = defaultopt32

# Include submodule files
include("radon.jl")
include("rjasanow.jl")

end # module

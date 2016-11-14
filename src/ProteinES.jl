module ProteinES

import Base: cos, sign, seek
import Base.LinAlg.BLAS: gemv!, axpy!
using Distances: euclidean

export
    # types.jl
    Element, Triangle, Tetrahedron, Charge, Option, SingleLayer, DoubleLayer, PotentialType, QuadraturePoints, QuadPts2D, QuadPts2D,
    # quad.jl
    quadraturepoints,
    # util.jl
    props!, meshunion,
    # this file
    defaultopt, ε0

# Global constants
const ε0 = 1 / (4π * 1e-7 * 299792458^2) # vacuum permittivity [F/m]

include("base/types.jl")
include("base/quad.jl")
include("base/util.jl")
include("base/molpot.jl")

# Default options
const defaultopt64 = Option(2., 78., 1.8, 20.)
const defaultopt32 = Option(2f0, 78f0, 1.8f0, 20f0)
defaultopt(::Type{Float64}) = defaultopt64
defaultopt(::Type{Float32}) = defaultopt32

# Include submodule files
include("IO.jl")
include("Radon.jl")
include("Rjasanow.jl")
include("BEM.jl")
#include("FEM.jl")

end # module

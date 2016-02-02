module ProteinES

import Base: cos, sign, seek
import Base.LinAlg.BLAS: gemv!, axpy!

export
    # types.jl
    Element, Triangle, Tetrahedron, Charge, Option, SingleLayer, DoubleLayer, PotentialType,
    # util.jl
    props!,
    # this file
    defaultopt, ε0

# Global constants
const ε0 = 1/ (4π * 1e-7 * 299792458^2)

include("base/types.jl")
include("base/util.jl")

# Default options
const defaultopt64 = Option(2., 78., 1.8, 20.)
const defaultopt32 = Option(2f0, 78f0, 1.8f0, 20f0)
defaultopt(::Type{Float64}) = defaultopt64
defaultopt(::Type{Float32}) = defaultopt32

# Include submodule files
include("IO.jl")
include("Radon.jl")
include("Rjasanow.jl")
include("Local.jl")
include("Nonlocal.jl")

end # module

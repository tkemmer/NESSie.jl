module ProteinES

import Base: cos, sign, seek
import Base.LinAlg.BLAS: gemv!, gemv, axpy!, gemm
using Distances: euclidean

include("base/model.jl")
export Element, Triangle, Tetrahedron, Charge, Model, SurfaceModel, VolumeModel

include("base/quadrature.jl")
export QuadraturePoints, QuadPts2D, QuadPts3D, quadraturepoints

include("base/constants.jl")
export ε0, σ, Option, defaultopt, yukawa

include("base/util.jl")
export props!, meshunion

include("base/potentials.jl")
export PotentialType, SingleLayer, DoubleLayer

# Submodules
include("IO.jl")
include("Radon.jl")
include("Rjasanow.jl")
include("BEM.jl")
#include("FEM.jl")

end # module

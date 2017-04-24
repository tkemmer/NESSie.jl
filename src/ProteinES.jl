module ProteinES

import Base: cos, sign, seek
import Base.LinAlg.BLAS: gemv!, gemv, axpy!, gemm
using Distances: euclidean

include("base/model.jl")
export Element, Triangle, Tetrahedron, Charge, Model, SurfaceModel, VolumeModel

include("base/quadrature.jl")
export QuadraturePoints, QuadPts2D, QuadPts3D, quadraturepoints

include("base/constants.jl")
export ε0, Option, defaultopt

include("base/util.jl")
export props!, meshunion, obspoints_line, obspoints_plane

include("base/potentials.jl")
export PotentialType, SingleLayer, DoubleLayer, LocalityType, LocalES, NonlocalES, φmol,
    ∂ₙφmol, ∇φmol

include("base/born.jl")
export BornIon, bornion, φΩ, φΣ

# Submodules
include("IO.jl")
include("Radon.jl")
include("Rjasanow.jl")
include("BEM.jl")
#include("FEM.jl")

end # module

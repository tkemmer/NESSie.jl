__precompile__()

module NESSie

import Base: cos, sign, seek
import LinearAlgebra.BLAS: gemv!, gemv, axpy!, gemm

using Distances: euclidean
using LinearAlgebra: ⋅, ×, norm, rmul!

include("base/constants.jl")
export ε0, Option, defaultopt

include("base/model.jl")
export Element, SurfaceElement, VolumeElement, Triangle, Tetrahedron, Charge, Model

include("base/quadrature.jl")
export QuadraturePoints, QuadPts2D, QuadPts3D, quadraturepoints

include("base/util.jl")
export meshunion, obspoints_line, obspoints_plane

include("base/potentials.jl")
export PotentialType, SingleLayer, DoubleLayer, LocalityType, LocalES, NonlocalES, φmol,
    ∂ₙφmol, ∇φmol

# Submodules
include("Format.jl")
export Format

include("Radon.jl")
export Radon

include("Rjasanow.jl")
export Rjasanow

include("BEM.jl")
export BEM

include("TestModel.jl")
export TestModel

end # module

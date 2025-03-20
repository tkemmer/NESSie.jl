__precompile__()

module NESSie

using AutoHashEquals
using Distances: euclidean
using Gmsh
using LinearAlgebra: ⋅, ×, BLAS, norm, rmul!

import FileIO, MeshIO
import GeometryBasics

include("base/constants.jl")
export ε0, Option, defaultopt

include("base/model.jl")
export Element, SurfaceElement, VolumeElement, Triangle, Tetrahedron, Charge, Model

include("base/quadrature.jl")
export QuadraturePoints, QuadPts2D, QuadPts3D, ElementQuad, TriangleQuad, quadraturepoints

include("base/util.jl")
export guess_domain, meshunion, obspoints_plane

include("base/potentials.jl")
export PotentialType, SingleLayer, DoubleLayer, LocalityType, LocalES, NonlocalES, φmol,
    ∂ₙφmol, ∇φmol, espotential, molpotential, rfenergy, rfpotential

include("base/deprecation.jl")

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

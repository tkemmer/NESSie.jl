#=
    Quadrature points
=#
abstract QuadraturePoints{T <: AbstractFloat}

#=
    Quadrature points and weights for triangles.

    Ref:
    [1] J. Radon. Zur mechanischen Kubatur. Monatsh. für Math. 52(4): 286-300, 1948.
=#
immutable QuadPts2D{T} <: QuadraturePoints{T}
    num::Int            # number of points
    x::Vector{T}        # x values
    y::Vector{T}        # y values
    weight::Vector{T}   # weights
end

const r15 = √15

const triquadpts64 = QuadPts2D(7,
    [1/3, (6+r15)/21, (9-2r15)/21, (6+r15)/21, (6-r15)/21, (9+2r15)/21, (6-r15)/21],
    [1/3, (6+r15)/21, (9-2r15)/21, (6+r15)/21, (6-r15)/21, (9+2r15)/21, (6-r15)/21],
    [9/80, (155+r15)/2400, (155+r15)/2400, (155+r15)/2400, (155-r15)/2400, (155-r15)/2400, (155-r15)/2400]
)

const triquadpts32 = QuadPts2D(7,
    Float32[1/3, (6+r15)/21, (9-2r15)/21, (6+r15)/21, (6-r15)/21, (9+2r15)/21, (6-r15)/21],
    Float32[1/3, (6+r15)/21, (9-2r15)/21, (6+r15)/21, (6-r15)/21, (9+2r15)/21, (6-r15)/21],
    Float32[9/80, (155+r15)/2400, (155+r15)/2400, (155+r15)/2400, (155-r15)/2400, (155-r15)/2400, (155-r15)/2400]
)

quadraturepoints(::Type{Triangle}, ::Type{Float64}) = triquadpts64
quadraturepoints(::Type{Triangle}, ::Type{Float32}) = triquadpts32

#=
    Quadrature points and weights for tetrahedra.

    Ref:
    [1] P Keast, Moderate degree tetrahedral quadrature formulas, CMAME 55: 339-348 (1986)
    TODO check
=#
immutable QuadPts3D{T} <: QuadraturePoints{T}
    num::Int            # number of points
    x::Vector{T}        # x values
    y::Vector{T}        # y values
    z::Vector{T}        # z values
    weight::Vector{T}   # weights
end

const tetraquadpts64 = QuadPts3D(5,
    [.25, .5, 1/6, 1/6, 1/6],
    [.25, 1/6, 1/6, 1/6, .5],
    [.25, 1/6, 1/6, .5, 1/6],
    [-.8, .45, .45, .45, .45]/6
)

const tetraquadpts32 = QuadPts3D(5,
    Float32[.25, .5, 1/6, 1/6, 1/6],
    Float32[.25, 1/6, 1/6, 1/6, .5],
    Float32[.25, 1/6, 1/6, .5, 1/6],
    Float32[-.8, .45, .45, .45, .45]/6f0
)

quadraturepoints(::Type{Tetrahedron}, ::Type{Float64}) = tetraquadpts64
quadraturepoints(::Type{Tetrahedron}, ::Type{Float32}) = tetraquadpts32

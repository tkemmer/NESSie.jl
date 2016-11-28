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
for T in [:Float64, :Float32]
    varname = Symbol("triquadpts_", T)
    @eval begin
        const $(varname) = QuadPts2D(7,
            $(T)[1/3, (6+r15)/21, (9-2r15)/21, (6+r15)/21, (6-r15)/21, (9+2r15)/21, (6-r15)/21],
            $(T)[1/3, (6+r15)/21, (9-2r15)/21, (6+r15)/21, (6-r15)/21, (9+2r15)/21, (6-r15)/21],
            $(T)[9/80, (155+r15)/2400, (155+r15)/2400, (155+r15)/2400, (155-r15)/2400, (155-r15)/2400, (155-r15)/2400]
        )
        quadraturepoints(::Type{Triangle}, ::Type{$(T)}) = $(varname)
    end
end

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

for T in [:Float64, :Float32]
    varname = Symbol("tetraquadpts_", T)
    @eval begin
        const $(varname) = QuadPts3D(5,
            $(T)[.25, .5, 1/6, 1/6, 1/6],
            $(T)[.25, 1/6, 1/6, 1/6, .5],
            $(T)[.25, 1/6, 1/6, .5, 1/6],
            $(T)[-.8, .45, .45, .45, .45]/6
        )
        quadraturepoints(::Type{Tetrahedron}, ::Type{$(T)}) = $(varname)
    end
end

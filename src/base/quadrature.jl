# =========================================================================================
"""
    abstract type QuadraturePoints{T <: AbstractFloat} end

Abstract base type for all sets of quadrature points
"""
abstract type QuadraturePoints{T <: AbstractFloat} end


# =========================================================================================
"""
    struct QuadPts2D{T} <: QuadraturePoints{T}
        num   ::Int         # number of points
        x     ::Vector{T}   # x values
        y     ::Vector{T}   # y values
        weight::Vector{T}   # weights
    end

Quadrature points and weights for triangles
"""
struct QuadPts2D{T} <: QuadraturePoints{T}
    """Number of points"""
    num::Int
    """x values"""
    x::Vector{T}
    """y values"""
    y::Vector{T}
    """Weights"""
    weight::Vector{T}
end


# =========================================================================================
"""
    struct QuadPts3D{T} <: QuadraturePoints{T}
        num   ::Int         # number of points
        x     ::Vector{T}   # x values
        y     ::Vector{T}   # y values
        z     ::Vector{T}   # z values
        weight::Vector{T}   # weights
    end

Quadrature points and weights for tetrahedra
"""
struct QuadPts3D{T} <: QuadraturePoints{T}
    """Number of points"""
    num::Int
    """x values"""
    x::Vector{T}
    """y values"""
    y::Vector{T}
    """z values"""
    z::Vector{T}
    """Weights"""
    weight::Vector{T}
end


# =========================================================================================
@doc """
    quadraturepoints(::Type{Tetrahedron}, ::Type{Float64})
    quadraturepoints(::Type{Tetrahedron}, ::Type{Float32})
    quadraturepoints(::Type{Triangle},    ::Type{Float64})
    quadraturepoints(::Type{Triangle},    ::Type{Float32})

Generator function for quadrature points:
 * *Triangles*: 7 points per element [[Rad48]](@ref Bibliography)
 * *Tetrahedra*: 5 points per element [[Kea86]](@ref Bibliography)

## Return type
`QuadPts2D` or `QuadPts3D`
""" quadraturepoints
for T in [:Float64, :Float32]
    varname = Symbol("triquadpts_", T)
    @eval begin
        const $(varname) = QuadPts2D(7,
            $(T)[1/3, (6+√15)/21, (9-2*√15)/21, (6+√15)/21,
                (6-√15)/21, (9+2*√15)/21, (6-√15)/21],
            $(T)[1/3, (6+√15)/21, (9-2*√15)/21, (6+√15)/21,
                (6-√15)/21, (9+2*√15)/21, (6-√15)/21],
            $(T)[9/80, (155+√15)/2400, (155+√15)/2400, (155+√15)/2400,
                (155-√15)/2400, (155-√15)/2400, (155-√15)/2400]
        )
        quadraturepoints(::Type{Triangle}, ::Type{$(T)}) = $(varname)
    end
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

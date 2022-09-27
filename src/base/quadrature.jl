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
    quadraturepoints(::Type{Tetrahedron{Float64}})
    quadraturepoints(::Type{Tetrahedron{Float32}})
    quadraturepoints(::Type{Triangle{Float64}})
    quadraturepoints(::Type{Triangle{Float32}})

Generator function for quadrature points:
 * *Triangles*: 7 points per element [[Rad48]](@ref Bibliography)
 * *Tetrahedra*: 5 points per element [[Kea86]](@ref Bibliography)

# Return type
`QuadPts2D` or `QuadPts3D`
""" quadraturepoints
for T in [:Float64, :Float32]
    varname = Symbol("triquadpts_", T)
    @eval begin
        const $(varname) = QuadPts2D(7,
            $(T)[1/3, (6+√15)/21, (9-2*√15)/21, (6+√15)/21,
                (6-√15)/21, (9+2*√15)/21, (6-√15)/21],
            $(T)[1/3, (9-2*√15)/21, (6+√15)/21, (6+√15)/21,
                (9+2*√15)/21, (6-√15)/21, (6-√15)/21],
            $(T)[9/80, (155+√15)/2400, (155+√15)/2400, (155+√15)/2400,
                (155-√15)/2400, (155-√15)/2400, (155-√15)/2400]
        )
        quadraturepoints(::Type{Triangle{$(T)}}) = $(varname)
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
        @inline quadraturepoints(::Type{Tetrahedron{$(T)}}) = $(varname)
    end
end


# =========================================================================================
"""
    abstract type ElementQuad{T <: AbstractFloat} end

Abstract base type for quadrature points on specific elements
"""
abstract type ElementQuad{T <: AbstractFloat} end


# =========================================================================================
"""
    struct TriangleQuad{T} <: ElementQuad{T}
        elem   ::Triangle{T}   # given element
        qpts   ::Matrix{T}     # quadrature points on the element
        weights::Vector{T}     # weights for the quadrature points
    end

Quadrature points on a specific surface triangle, including weights.

# Special constructors
```julia
TriangleQuad{T}(elem::Triangle{T})
````
Computes quadrature points and weights for the given triangle.
"""
struct TriangleQuad{T} <: ElementQuad{T}
    """Given element"""
    elem   ::Triangle{T}
    """Quadrature points on the given element as 3xN matrix"""
    qpts   ::Matrix{T}
    """Weights for the quadrature points"""
    weights::Vector{T}
end

@inline function TriangleQuad(elem::Triangle{T}) where T
    qpts = quadraturepoints(Triangle{T})
    mat  = Matrix{T}(undef, 3, qpts.num)
    for j in 1:qpts.num
        mat[:, j] .= qpts.x[j] .* (elem.v2 .- elem.v1) .+ qpts.y[j] .* (elem.v3 .- elem.v1) .+ elem.v1
    end
    TriangleQuad{T}(elem, mat, qpts.weight)
end


# =========================================================================================
"""
    quadraturepoints{T}(elements::Vector{Triangle{T}})

Computes and returns quadrature points on all given elements, including weights.

# Return type
`Vector{TriangleQuad{T}}`
"""
function quadraturepoints(elements::Vector{Triangle{T}}) where T
    qpts = quadraturepoints(Triangle{T})
    ret  = Vector{TriangleQuad{T}}(undef, length(elements))

    for i in eachindex(elements)
        elem = elements[i]
        mat  = Matrix{T}(undef, 3, qpts.num)
        for j in 1:qpts.num
            mat[:, j] .= qpts.x[j] .* (elem.v2 .- elem.v1) .+ qpts.y[j] .* (elem.v3 .- elem.v1) .+ elem.v1
        end
        ret[i] = TriangleQuad{T}(elem, mat, qpts.weight)
    end
    ret
end

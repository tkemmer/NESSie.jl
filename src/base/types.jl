#=
    Elements
=#
abstract Element{T <: AbstractFloat}

#=
    Representation of a single surface triangle
=#
type Triangle{T} <: Element{T}
    v1::Vector{T}       # position of the first node
    v2::Vector{T}       # position of the second node
    v3::Vector{T}       # position of the third node
    center::Vector{T}   # centroid of the triangle
    normal::Vector{T}   # normal vector of the triangle
    area::T             # area of the triangle
    distorig::T         # distance to the origin
end
Triangle{T}(v1::Vector{T}, v2::Vector{T}, v3::Vector{T}) = Triangle(v1, v2, v3, T[], T[], zero(T), zero(T))

#=
    Representation of a tetrahedron
=#
type Tetrahedron{T} <: Element{T}
    v1::Vector{T}       # position of the first node
    v2::Vector{T}       # position of the second node
    v3::Vector{T}       # position of the third node
    v4::Vector{T}       # position of the fourth node
    domain::Symbol      # element domain (solvent :Σ, solute :Ω, or :none)
end
Tetrahedron{T}(v1::Vector{T}, v2::Vector{T}, v3::Vector{T}, v4::Vector{T}) = Tetrahedron(v1, v2, v3, v4, :none)

#=
    Representation of a single point charge
=#
type Charge{T <: AbstractFloat}
    pos::Vector{T}  # position of the charge
    val::T          # charge value
end
Charge{T <: AbstractFloat}(args::T...) = Charge{T}([args[1:end-1]...], args[end])
Charge(T::DataType, args...) = Charge{T}(convert(Vector{T}, [args[1:end-1]...]), convert(T, args[end]))

#=
    Constants
    TODO move yukawa to a more suitable place
=#
immutable Option{T <: AbstractFloat}
    εΩ::T       # dielectric constant of the solute
    εΣ::T       # dielectric constant of the solvent
    ε∞::T       # large-scale (bulk) solvent response
    λ::T        # scale
    yukawa::T   # exponent for fundamental solution of yukawa operator -1/Λ = -1/(λ√(ε∞/εΣ))
end
Option{T <: AbstractFloat}(εΩ::T, εΣ::T, ε∞::T, λ::T) = Option(εΩ, εΣ, ε∞, λ, √(εΣ/ε∞)/λ)

#=
    Enum-like representation of potential types
=#
abstract PotentialType
type SingleLayer <: PotentialType end
type DoubleLayer <: PotentialType end

#=
    Quadrature points
=#
abstract QuadraturePoints{T <: AbstractFloat}

#=
    Quadrature points in 2D
=#
immutable QuadPts2D{T} <: QuadraturePoints{T}
    num::Int            # number of points
    x::Vector{T}        # x values
    y::Vector{T}        # y values
    weight::Vector{T}   # weights
end

#=
    Quadrature points in 3D
=#
immutable QuadPts3D{T} <: QuadraturePoints{T}
    num::Int            # number of points
    x::Vector{T}        # x values
    y::Vector{T}        # y values
    z::Vector{T}        # z values
    weight::Vector{T}   # weights
end

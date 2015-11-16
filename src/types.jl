abstract Element{T <: AbstractFloat}

#=
    Representation of a single surface triangle.
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
    Representation of a single charge.
=#
type Charge{T <: AbstractFloat}
    pos::Vector{T}  # position of the charge
    val::T          # charge value
end
Charge{T <: AbstractFloat}(args::T...) = Charge{T}([args[1:end-1]...], args[end])
Charge(dtype::DataType, args...) = Charge{dtype}(convert(Vector{dtype}, [args[1:end-1]...]), convert(dtype, args[end]))

#=
    Constants.
=#
immutable Option{T <: AbstractFloat}
    εΩ::T       # dielectric constant of the solute
    εΣ::T       # dielectric constant of the solvent
    ε∞::T       # large-scale (bulk) solvent response
    λ::T        # scale
    yukawa::T   # exponent for fundamental solution of yukawa operator -1/Λ = -1/(λ√(ε∞/εΣ))
end
Option{T <: AbstractFloat}(εΩ::T, εΣ::T, ε∞::T, λ::T) = Option(εΩ, εΣ, ε∞, λ, -√(εΣ/ε∞)/λ)

#=
    Enum-like representation of potential types.
=#
abstract PotentialType
type SingleLayer <: PotentialType end
type DoubleLayer <: PotentialType end

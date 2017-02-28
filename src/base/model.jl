#=
    Elements
=#
abstract type Element{T <: AbstractFloat} end

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
    Models
=#
abstract type Model{T <: AbstractFloat} end

#=
    Surface model
    Protein immersed in a structured solvent, represented by the protein surface and its partial point charges.
=#
type SurfaceModel{T} <: Model{T}
    elements::Vector{Triangle{T}}
    charges::Vector{Charge{T}}
end

#=
    Volume model
    Protein immersed in a structured solvent; protein and a sphere of surrounding space represented as a collection of
    tetrahedra and the protein's partial point charges.
=#
type VolumeModel{T} <: Model{T}
    elements::Vector{Tetrahedron{T}}
    charges::Vector{Charge{T}}
end

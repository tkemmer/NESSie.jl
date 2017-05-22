# =========================================================================================
"""
    abstract type Element{T <: AbstractFloat} end

Abstract base type for all elements.
"""
abstract type Element{T <: AbstractFloat} end


# =========================================================================================
"""
    type Triangle{T} <: Element{T}
        v2      ::Vector{T}   # position of the second node
        v1      ::Vector{T}   # position of the first node
        v3      ::Vector{T}   # position of the third node
        center  ::Vector{T}   # centroid of the triangle
        normal  ::Vector{T}   # normal vector of the triangle
        area    ::T           # area of the triangle
        distorig::T           # distance to the origin
    end

Representation of a single surface triangle.

# Special constructors
```julia
Triangle{T}(
    v1  ::Vector{T},
    v2  ::Vector{T},
    v3  ::Vector{T},
    init::Bool=true
)
```
Most commonly used constructor variant for creating a triangle by only specifying its nodes.
If `init` is set to `true`, the remaining member variables will automatically be computed
via [`props!`](@ref ProteinES.props!) to fully initialize the object.

!!! warning
    Most ProteinES.jl functions assume that given triangles are fully initialized. Using
    these functions with partly initialized triangles can lead to undefined behavior!
"""
type Triangle{T} <: Element{T}
    "position of the first node"
    v1::Vector{T}
    "position of the second node"
    v2::Vector{T}
    "position of the third node"
    v3::Vector{T}
    "centroid of the triangle"
    center::Vector{T}
    "normal vector of the triangle"
    normal::Vector{T}
    "area of the triangle"
    area::T
    "distance to the origin"
    distorig::T
end
Triangle{T}(v1::Vector{T}, v2::Vector{T}, v3::Vector{T}, init::Bool=true) = begin
    r = Triangle(v1, v2, v3, T[], T[], zero(T), zero(T))
    init && props!(r)
    r
end


# =========================================================================================
"""
    type Tetrahedron{T} <: Element{T}
        v1::Vector{T}       # position of the first node
        v2::Vector{T}       # position of the second node
        v3::Vector{T}       # position of the third node
        v4::Vector{T}       # position of the fourth node
        domain::Symbol      # element domain (solvent :Σ, solute :Ω, or :none)
    end

Representation of a single tetrahedron.

# Special constructors
```julia
Tetrahedron{T}(
    v1::Vector{T},
    v2::Vector{T},
    v3::Vector{T},
    v4::Vector{T}
)
```
Sets domain to `:none`.
"""
type Tetrahedron{T} <: Element{T}
    "position of the first node"
    v1::Vector{T}
    "position of the second node"
    v2::Vector{T}
    "position of the third node"
    v3::Vector{T}
    "position of the fourth node"
    v4::Vector{T}
    "element domain (solvent :Σ, solute :Ω, or :none)"
    domain::Symbol
end
Tetrahedron{T}(v1::Vector{T}, v2::Vector{T}, v3::Vector{T}, v4::Vector{T}) =
    Tetrahedron(v1, v2, v3, v4, :none)


# =========================================================================================
"""
    type Charge{T <: AbstractFloat}
        pos::Vector{T}  # position of the charge
        val::T          # charge value
    end

Representation of a single point charge.

# Special constructors
```julia
Charge{T}(
    posx::T,
    posy::T,
    posz::T,
    val ::T
)
```
Constructor variant with flat argument list for `pos`.
"""
type Charge{T <: AbstractFloat}
    "position of the charge"
    pos::Vector{T}
    "charge value"
    val::T
end
Charge{T}(posx::T, posy::T, posz::T, val::T) = Charge{T}([posx, posy, posz], val)


# =========================================================================================
"""
    abstract type Model{T <: AbstractFloat} end

Abstract base type for all models.
"""
abstract type Model{T <: AbstractFloat} end


# =========================================================================================
"""
    type SurfaceModel{T} <: Model{T}
        nodes   ::Vector{Vector{T}    = Vector{T}[]
        elements::Vector{Triangle{T}} = Triangle{T}[]
        charges ::Vector{Charge{T}}   = Charge{T}[]
    end

Surface model; Typically a protein immersed in a structured solvent, represented by the
protein surface and its partial point charges.
"""
type SurfaceModel{T} <: Model{T}
    nodes   ::Vector{Vector{T}}
    elements::Vector{Triangle{T}}
    charges ::Vector{Charge{T}}

    SurfaceModel{T}(
        nodes   ::Vector{Vector{T}}   = Vector{T}[],
        elements::Vector{Triangle{T}} = Triangle{T}[],
        charges ::Vector{Charge{T}}   = Charge{T}[]
    ) where T = new(nodes, elements, charges)
end

SurfaceModel(
    nodes   ::Vector{Vector{T}}   = Vector{T}[],
    elements::Vector{Triangle{T}} = Triangle{T}[],
    charges ::Vector{Charge{T}}   = Charge{T}[],
) where T = SurfaceModel{T}(nodes, elements, charges)


# =========================================================================================
"""
    type VolumeModel{T} <: Model{T}
        nodes   ::Vector{Vector{T}}      = Vector{T}[]
        elements::Vector{Tetrahedron{T}} = Tetrahedron{T}[]
        charges ::Vector{Charge{T}}      = Charge{T}[]
    end

Volume model; Typically a protein immersed in a structured solvent. The protein and a sphere
of surrounding space are represented as a collection of tetrahedra and the protein's partial
point charges.
"""
type VolumeModel{T} <: Model{T}
    nodes   ::Vector{Vector{T}}
    elements::Vector{Tetrahedron{T}}
    charges ::Vector{Charge{T}}
    params  ::Option{T}

    VolumeModel{T}(
        nodes   ::Vector{Vector{T}}      = Vector{T}[],
        elements::Vector{Tetrahedron{T}} = Tetrahedron{T}[],
        charges ::Vector{Charge{T}}      = Charge{T}[]
    ) where T = new(nodes, elements, charges)
end

VolumeModel(
    nodes   ::Vector{Vector{T}}      = Vector{T}[],
    elements::Vector{Tetrahedron{T}} = Tetrahedron{T}[],
    charges ::Vector{Charge{T}}      = Charge{T}[]
) where T = VolumeModel{T}(nodes, elements, charges)

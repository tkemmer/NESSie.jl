# =========================================================================================
"""
    abstract type Element{T <: AbstractFloat} end
    abstract type SurfaceElement{T} <: Element{T} end
    abstract type VolumeElement{T}  <: Element{T} end

Abstract base types for all elements.
"""
abstract type Element{T <: AbstractFloat} end
abstract type SurfaceElement{T} <: Element{T} end
abstract type VolumeElement{T}  <: Element{T} end


# =========================================================================================
"""
    struct Triangle{T} <: SurfaceElement{T}
        v1      ::Vector{T}   # position of the first node
        v2      ::Vector{T}   # position of the second node
        v3      ::Vector{T}   # position of the third node
        center  ::Vector{T}   # centroid of the triangle
        normal  ::Vector{T}   # normal vector of the triangle
        area    ::T           # area of the triangle
        distorig::T           # distance to the origin
    end

Representation of a single surface triangle.

# Special constructors
```julia
Triangle(
    v1  ::Vector{T},
    v2  ::Vector{T},
    v3  ::Vector{T}
)
```
Most commonly used constructor variant for creating a triangle by only specifying its nodes.
The remaining member variables will automatically be computed via
[`props`](@ref NESSie.props).
"""
@auto_hash_equals struct Triangle{T} <: SurfaceElement{T}
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

@inline Triangle(
    v1::Vector{T},
    v2::Vector{T},
    v3::Vector{T}
) where T = props(Triangle{T}(v1, v2, v3, T[], T[], zero(T), zero(T)))

@inline function Base.show(io::IO, ::MIME"text/plain", elem::Triangle)
    show(io, elem)
end

@inline function Base.show(io::IO, elem::Triangle)
    print(io,
        "$(typeof(elem))",
        "(v1 = $(repr(elem.v1))",
        ", v2 = $(repr(elem.v2))",
        ", v3 = $(repr(elem.v3)))",
    )
end


# =========================================================================================
"""
    struct Tetrahedron{T} <: VolumeElement{T}
        v1::Vector{T}       # position of the first node
        v2::Vector{T}       # position of the second node
        v3::Vector{T}       # position of the third node
        v4::Vector{T}       # position of the fourth node
        domain::Symbol      # element domain (solvent :Σ, solute :Ω, or :none)
    end

Representation of a single tetrahedron.

# Special constructors
```julia
Tetrahedron(
    v1::Vector{T},
    v2::Vector{T},
    v3::Vector{T},
    v4::Vector{T}
)
```
Sets domain to `:none`.
"""
@auto_hash_equals struct Tetrahedron{T} <: VolumeElement{T}
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

@inline Tetrahedron(
    v1::Vector{T},
    v2::Vector{T},
    v3::Vector{T},
    v4::Vector{T}
) where T = Tetrahedron{T}(v1, v2, v3, v4, :none)

@inline function Base.show(io::IO, ::MIME"text/plain", elem::Tetrahedron)
    show(io, elem)
end

@inline function Base.show(io::IO, elem::Tetrahedron)
    print(io,
        "$(typeof(elem))",
        "(v1 = $(repr(elem.v1))",
        ", v2 = $(repr(elem.v2))",
        ", v3 = $(repr(elem.v3))",
        ", v4 = $(repr(elem.v4)))",
    )
end


# =========================================================================================
"""
    struct Charge{T <: AbstractFloat}
        pos::Vector{T}  # position of the charge
        val::T          # charge value
    end

Representation of a single point charge.

# Special constructors
```julia
Charge(
    posx::T,
    posy::T,
    posz::T,
    val ::T
)
```
Constructor variant with flat argument list for `pos`.
"""
@auto_hash_equals struct Charge{T <: AbstractFloat}
    "position of the charge"
    pos::Vector{T}
    "charge value"
    val::T
end

@inline Charge(
    posx::T,
    posy::T,
    posz::T,
    val ::T
) where T = Charge{T}([posx, posy, posz], val)

@inline function Base.show(io::IO, ::MIME"text/plain", charge::Charge)
    show(io, charge)
end

@inline function Base.show(io::IO, charge::Charge)
    print(io,
        "$(typeof(charge))",
        "(pos = $(repr(charge.pos))",
        ", val = $(charge.val))"
    )
end


# =========================================================================================
"""
    mutable struct Model{T, E <: Element{T}}
        nodes   ::Vector{Vector{T}  = Vector{T}[]    # mesh nodes
        elements::Vector{E}         = E[]            # mesh elements
        charges ::Vector{Charge{T}} = Charge{T}[]    # point charges in the molecule
        params  ::Option{T}         = defaultopt(T)  # system constants
    end

System model representing a biomelecule in solvation, including a collection of point
charges in the molecule and a set of system constants. The system can either be represented
as a surface model (e.g., a collection of molecule surface triangles) or as a volume model
(e.g., a collection of tetrahedra for the molecule and its surrounding space).
"""
@auto_hash_equals mutable struct Model{T, E <: Element{T}}
    """mesh nodes"""
    nodes   ::Vector{Vector{T}}
    """mesh elements"""
    elements::Vector{E}
    """point charges in the molecule"""
    charges ::Vector{Charge{T}}
    """system constants"""
    params  ::Option{T}

    Model{T, E}(
        nodes   ::Vector{Vector{T}} = Vector{T}[],
        elements::Vector{E}         = E[],
        charges ::Vector{Charge{T}} = Charge{T}[],
        params  ::Option{T}         = defaultopt(T)
    ) where {T, E <: Element{T}}    = new(nodes, elements, charges, params)
end

@inline Model(
    nodes   ::Vector{Vector{T}},
    elements::Vector{E},
    charges ::Vector{Charge{T}} = Charge{T}[],
    params  ::Option{T}         = defaultopt(T)
) where {T, E <: Element{T}}    = Model{T, E}(nodes, elements, charges, params)

@inline function Base.show(io::IO, ::MIME"text/plain", model::Model)
    show(io, model)
end

@inline function Base.show(io::IO, model::Model)
    print(io,
        "$(typeof(model))",
        "(nodes = ", length(model.nodes),
        ", elements = ", length(model.elements),
        ", charges = ", length(model.charges), ")"
    )
end

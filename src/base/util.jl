"""
    props!{T}(
        elem::Triangle{T}
    )

Computes the given triangle's properties, that is, centroid, normal, distance to origin,
and area.

# Return type
`Void`
"""
function props!{T}(elem::Triangle{T})
    # reject degenerate triangles
    @assert !isdegenerate(elem) "Degenerate triangle $(elem)"

    # compute centroid
    elem.center = 3 \ (elem.v1 + elem.v2 + elem.v3)

    # compute normal
    elem.normal = (elem.v2 - elem.v1) × (elem.v3 - elem.v1)
    vnorm = vecnorm(elem.normal)
    elem.normal /= vnorm

    # compute distance to origin
    elem.distorig = elem.normal ⋅ elem.v1

    # compute area
    elem.area = 2 \ vnorm
    nothing
end


#TODO VolumeModel
"""
    meshunion{T}(
        nodesΩ   ::Vector{Vector{T}},
        elementsΩ::Vector{Tetrahedron{T}},
        nodesΣ   ::Vector{Vector{T}},
        elementsΣ::Vector{Tetrahedron{T}}
    )

Merges two lists of nodes and elements, e.g., the volume meshes of a protein and the
solvent. Duplicate nodes (e.g., the nodes on the protein surface) are merged into a single
node, duplicate elements (if any) are retained.

# Return type
`(Vector{Vector{T}}, Vector{Tetrahedron{T}})`

!!! note
    This function assumes that there are no duplicates within either of the node lists!
"""
function meshunion{T}(
        nodesΩ::Vector{Vector{T}},
        elementsΩ::Vector{Tetrahedron{T}},
        nodesΣ::Vector{Vector{T}},
        elementsΣ::Vector{Tetrahedron{T}}
    )
    # find nodes that are to be replaced (tbr)
    obsolete = nodesΣ ∩ nodesΩ
    tbr = Dict{Vector{T}, Int}(zip(obsolete, indexin(obsolete, nodesΩ)))

    elements = copy(elementsΣ)
    for elem in elementsΣ
        haskey(tbr, elem.v1) && (elem.v1 = nodesΩ[tbr[elem.v1]])
        haskey(tbr, elem.v2) && (elem.v2 = nodesΩ[tbr[elem.v2]])
        haskey(tbr, elem.v3) && (elem.v3 = nodesΩ[tbr[elem.v3]])
        haskey(tbr, elem.v4) && (elem.v4 = nodesΩ[tbr[elem.v4]])
    end
    append!(elements, elementsΩ)

    collect(Set(nodesΣ) ∪ Set(nodesΩ)), elements
end


"""
    unpack{T}(data::Vector{Vector{T}})

Unpacks the given vector of vectors into a single vector.

# Return type
`Vector{T}`

# Example

```jldoctest
julia> unpack([[1, 2], [3]])
3-element Array{Int64,1}:
 1
 2
 3
```
"""
function unpack{T}(data::Vector{Vector{T}})
    isempty(data) ?  T[] : T[x for y in data for x in y]
end


"""
    vertexnormals{T}(model::SurfaceModel{T})

Returns a vector containing the normal vectors of the given surface model's nodes.

# Return type
`Vector{Vector{T}}`
"""
function vertexnormals{T}(model::SurfaceModel{T})
    revidx = reverseindex(model.nodes)
    normals = Vector{T}[zeros(T, 3) for _ in 1:length(model.nodes)]
    count = zeros(T, length(model.nodes))
    @inbounds for elem in model.elements, node in (elem.v1, elem.v2, elem.v3)
        idx = revidx[object_id(node)]
        count[idx] += 1
        for i in 1:3
            normals[idx][i] += (elem.normal[i]-normals[idx][i]) / count[idx]
        end
    end
    normals
end


"""
    eye!{T}(
        m::Union{DenseArray{T,2}, SubArray{T,2}},
        α::Number=one(T)
    )

Initializes the given matrix `m` with `αI`, with `I` being an identity matrix with the same
dimensions as `m`.

# Return type
`Void`

# Example
```jldoctest
julia> m = 2 * ones(2, 2)
2×2 Array{Float64,2}:
 2.0  2.0
 2.0  2.0

julia> eye!(m); m
2×2 Array{Float64,2}:
 1.0  0.0
 0.0  1.0

julia> eye!(m, 2); m
2×2 Array{Float64,2}:
 2.0  0.0
 0.0  2.0
```
"""
function eye!{T}(m::Union{DenseArray{T,2}, SubArray{T,2}}, α::Number=one(T))
    fill!(m, zero(T))
    pluseye!(m, α)
end


"""
    pluseye!{T}(
        m::Union{DenseArray{T,2}, SubArray{T,2}},
        α::Number=one(T)
    )

Adds `α` to all diagonal elements of matrix `m`.

# Return type
`Void`

# Example
```jldoctest
julia> m = 2 * ones(2, 2)
2×2 Array{Float64,2}:
 2.0  2.0
 2.0  2.0

julia> pluseye!(m); m
2×2 Array{Float64,2}:
 3.0  2.0
 2.0  3.0

julia> pluseye!(m, 2); m
2×2 Array{Float64,2}:
 5.0  2.0
 2.0  5.0
```
"""
function pluseye!{T}(m::Union{DenseArray{T,2}, SubArray{T,2}}, α::Number=one(T))
    α = convert(T, α)
    @inbounds for i in 1:min(size(m)...)
        m[i, i] += α
    end
    nothing
end


"""
    isdegenerate{T}(elem::Triangle{T})

Tests whether the given triangle is degenerate.

# Return type
`Bool`
"""
function isdegenerate{T}(elem::Triangle{T})
    @assert length(elem.v1) == length(elem.v2) == length(elem.v3) == 3
    u1 = elem.v2 - elem.v1
    u2 = elem.v3 - elem.v1
    cosine = u1 ⋅ u2 / vecnorm(u1) / vecnorm(u2)
    vecnorm(elem.v1 - elem.v2) < 1e-10 || vecnorm(elem.v1 - elem.v3) < 1e-10 ||
    vecnorm(elem.v2 - elem.v3) < 1e-10 || 1 - abs(cosine) <= 1e-10
end


"""
    seek(
        fh         ::IOStream,
        prefix     ::String,
        skiptheline::Bool = true
    )

Fast-forwards an IOStream to the next line starting with the given `prefix`. In case there
is no such line, the stream handle will be set to EOF.

# Arguments
 * `skiptheline`    If true, said line will also be skipped

# Return type
`Void`
"""
function seek(fh::IOStream, prefix::String, skiptheline::Bool=true)
    m = -1
    while !eof(fh)
        skiptheline || (m = position(fh))
        startswith(readline(fh), prefix) && break
    end
    eof(fh) || skiptheline || seek(fh, m)
    nothing
end


"""
    cos{T}(
        u    ::Vector{T},
        v    ::Vector{T},
        unorm::T=vecnorm(u),
        vnorm::T=vecnorm(v)
    )

Computes the cosine of the angle between the given vectors `u` and `v` with lengths `unorm`
and `vnorm`, respectively.

# Return type
`T`
"""
function cos{T}(u::Vector{T}, v::Vector{T}, unorm::T=vecnorm(u), vnorm::T=vecnorm(v))
    u ⋅ v / (unorm * vnorm)
end


"""
    cathetus{T}(hyp::T, cosθ::T)

Computes the cathetus ``c₁`` of a triangle given the hypotenuse ``h`` and the cosine of the
exterior angle ``θ`` between the hypotenuse and the other cathetus ``c₂``.

```math
c₂ = h ⋅ \\cos(θ) \\\\
h² = c₁² + c₂² \\\\
⇔ c₁ = \\sqrt{h² ⋅ (1 - \\cos²(θ))}
```

# Return type
`T`
"""
function cathetus{T}(hyp::T, cosθ::T)
    √(hyp^2 * (1 - cosθ^2))
end


"""
    sign{T}(
        u::Vector{T},
        v::Vector{T},
        n::Vector{T}
    )

Determines whether the normal vector of the plane specified by the vectors `u` and `v` has
the same orientation as the given normal vector `n`. Returns ``1`` if both normals have the
same orientation, ``0`` if at least one of the vectors is zero, and ``-1`` otherwise.

# Return type
`T`
"""
function sign{T}(u::Vector{T}, v::Vector{T}, n::Vector{T})
    # Devectorized version of sign((u1 × u2) ⋅ normal)
    sign(
        (u[2]*v[3] - u[3]*v[2]) * n[1] +
        (u[3]*v[1] - u[1]*v[3]) * n[2] +
        (u[1]*v[2] - u[2]*v[1]) * n[3]
    )
end


"""
    distance{T}(
        q   ::Vector{T},
        elem::Triangle{T}
    )

Calculates the (positive or negative) distance from the given point `q` to the plane the
given triangle `elem` is located in.

# Return type
`T`
"""
function distance{T}(q::Vector{T}, elem::Triangle{T})
    q ⋅ elem.normal - elem.distorig
end


"""
    ddot{T}(
        u::Vector{T},
        v::Vector{T},
        n::Vector{T}
    )

Devectorized computation of `(u-v)⋅n`.

# Return type
`T`
"""
function ddot{T}(u::Vector{T}, v::Vector{T}, n::Vector{T})
    (u[1] - v[1]) * n[1] + (u[2] - v[2]) * n[2] + (u[3] - v[3]) * n[3]
end


"""
    reverseindex{T}(v::Vector{T})

Creates a reverse index for the given vector `v`, that is, a dictionary linking the object
IDs of the vector elements to the corresponding position in the vector.

# Return type
`Dict{UInt, UInt}`
"""
function reverseindex{T}(v::Vector{T})
    Dict{UInt, UInt}(object_id(e) => i for (i,e) in enumerate(v))
end


"""
    obspoints_line{T}(
        u::Vector{T},
        v::Vector{T},
        n::Int
    )

Generates `n` evenly distributed observation points along the line segment from `u` to `v`.

# Return type
`Function`

# Example
```julia
for ξ in obspoints_line([0, 0, 0], [1, 1, 1], 10)
    ...
end
```
"""
function obspoints_line{T}(u::Vector{T}, v::Vector{T}, n::Int)
    (u + T(i) * (v - u) for i in linspace(0, 1, n))
end


"""
    obspoints_plane{T}(
        a  ::Vector{T},
        b  ::Vector{T},
        c  ::Vector{T},
        nba::Int,
        nbc::Int
    )

Generates `nba` ⋅ `nbc` evenly distributed observation points on the parallelogram with the
sides ``\\overline{BA}`` and ``\\overline{BC}`` and `a`, `b`, `c` being the location vectors
to the points ``A``, ``B``, and ``C``, respectively.

# Arguments
 * `nba`   Number of observation points along ``\\overline{BA}``
 * `nbc`   Number of observation points along ``\\overline{BC}``

# Return type
`Function`

# Example
```julia
for Ξ in obspoints_plane(...)
    for ξ in Ξ
        ...
    end
end
```
"""
function obspoints_plane{T}(a::Vector{T}, b::Vector{T}, c::Vector{T}, nba::Int, nbc::Int)
    (obspoints_line(ξ, c + (ξ - b), nbc) for ξ in obspoints_line(b, a, nba))
end


# Convenience aliases
gemv!{T}(α::T, m::Union{DenseArray{T,2}, SubArray{T,2}}, v::Vector{T}, dest::Union{DenseArray{T,1}, SubArray{T,1}}) = gemv!(α, m, v, one(T), dest)
gemv!{T}(α::T, m::Union{DenseArray{T,2}, SubArray{T,2}}, v::Vector{T}, β::T, dest::Union{DenseArray{T,1}, SubArray{T,1}}) = gemv!('N', α, m, v, β, dest)
gemv{T}(α::T, m::Union{DenseArray{T,2}, SubArray{T,2}}, v::Vector{T}) = gemv('N', α, m, v)
gemm{T}(α::T, a::Union{DenseArray{T,2}, SubArray{T, 2}}, b::Union{DenseArray{T,2}, SubArray{T,2}}) = gemm('N', 'N', α, a, b)

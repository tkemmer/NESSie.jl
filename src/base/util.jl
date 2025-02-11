# =========================================================================================
"""
    props(
        elem::Triangle{T}
    )

Computes the given triangle's properties, that is, centroid, normal, distance to origin,
and area. Returns the completely initialized Triangle as a copy.

!!! warning
    The given triangle remains unchanged!

# Return type
`Triangle`
"""
function props(elem::Triangle)
    # reject degenerate triangles
    @assert !isdegenerate(elem) "Degenerate triangle $(elem)"

    # compute centroid
    center = (elem.v1 .+ elem.v2 .+ elem.v3) ./ 3

    # compute normal
    normal = (elem.v2 .- elem.v1) × (elem.v3 .- elem.v1)
    vnorm = _norm(normal)
    normal ./= vnorm

    # compute distance to origin
    distorig = _dot(normal, elem.v1)

    # compute area
    area = vnorm / 2

    Triangle(elem.v1, elem.v2, elem.v3, center, normal, area, distorig)
end


# =========================================================================================
"""
    meshunion(
        model1::Model{T, Tetrahedron{T}},
        model2::Model{T, Tetrahedron{T}}
    )

Merges two volume models, e.g., the models of a protein and the solvent. Duplicate nodes
(e.g., the nodes on the protein surface) are merged into a single node, duplicate elements
and charges (if any) are retained as well as the system constants of the first model.

# Return type
`Model{T, Tetrahedron{T}}`

!!! note
    This function assumes that there are no duplicates within either of the node lists!
"""
function meshunion(
        model1::Model{T, Tetrahedron{T}},
        model2::Model{T, Tetrahedron{T}}
    ) where T
    # find nodes that are to be replaced (tbr)
    obsolete = model1.nodes ∩ model2.nodes
    indices = something.(indexin(obsolete, model1.nodes)) # Vector{Union{Nothing, Int}} -> Vector{Int}
    tbr = Dict{Vector{T}, Int}(zip(obsolete, indices))

    elements = copy(model1.elements)
    for elem in model2.elements
        v1 = haskey(tbr, elem.v1) ? model1.nodes[tbr[elem.v1]] : elem.v1
        v2 = haskey(tbr, elem.v2) ? model1.nodes[tbr[elem.v2]] : elem.v2
        v3 = haskey(tbr, elem.v3) ? model1.nodes[tbr[elem.v3]] : elem.v3
        v4 = haskey(tbr, elem.v4) ? model1.nodes[tbr[elem.v4]] : elem.v4
        push!(elements, Tetrahedron(v1, v2, v3, v4, elem.domain))
    end

    Model(
        collect(Vector{T}, Set(model2.nodes) ∪ Set(model1.nodes)),
        elements,
        model1.charges ∪ model2.charges,
        model1.params
    )
end


# =========================================================================================
"""
    unpack(data::Vector{Vector{T}})

Unpacks the given vector of vectors into a single vector.

# Return type
`Vector{T}`

# Example

```jldoctest; setup = :(using NESSie: unpack)
julia> unpack([[1, 2], [3]])
3-element Vector{Int64}:
 1
 2
 3
```
"""
@inline function unpack(data::AbstractVector{Vector{T}}) where T
    collect(T, x for y in data for x in y)
end


# =========================================================================================
"""
    vertexnormals(model::Model{T, Triangle{T}})

Returns a vector containing the normal vectors of the given model's triangles.

# Return type
`Vector{Vector{T}}`
"""
function vertexnormals(model::Model{T, Triangle{T}}) where T
    revidx = _reverseindex(model.nodes)
    normals = Vector{T}[zeros(T, 3) for _ in 1:length(model.nodes)]
    count = zeros(T, length(model.nodes))
    @inbounds for elem in model.elements, node in (elem.v1, elem.v2, elem.v3)
        idx = revidx[node]
        count[idx] += 1
        normals[idx] .+= (elem.normal .- normals[idx]) ./ count[idx]
    end
    normals
end


# =========================================================================================
"""
    eye!(
        m::AbstractMatrix{T},
        α::Number=one(T)
    )

Initializes the given matrix `m` with `αI`, with `I` being an identity matrix with the same
dimensions as `m`.

# Return type
`Void`

# Example
```jldoctest; setup = :(using NESSie: eye!)
julia> m = 2 * ones(2, 2)
2×2 Matrix{Float64}:
 2.0  2.0
 2.0  2.0

julia> eye!(m); m
2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0

julia> eye!(m, 2); m
2×2 Matrix{Float64}:
 2.0  0.0
 0.0  2.0
```
"""
@inline function eye!(
        m::AbstractMatrix{T},
        α::Number=one(T)
    ) where T
    fill!(m, zero(T))
    pluseye!(m, α)
end


# =========================================================================================
"""
    pluseye!(
        m::AbstractMatrix{T},
        α::Number=one(T)
    )

Adds `α` to all diagonal elements of matrix `m`.

# Return type
`Void`

# Example
```jldoctest; setup = :(using NESSie: pluseye!)
julia> m = 2 * ones(2, 2)
2×2 Matrix{Float64}:
 2.0  2.0
 2.0  2.0

julia> pluseye!(m); m
2×2 Matrix{Float64}:
 3.0  2.0
 2.0  3.0

julia> pluseye!(m, 2); m
2×2 Matrix{Float64}:
 5.0  2.0
 2.0  5.0
```
"""
function pluseye!(
        m::AbstractMatrix{T},
        α::Number=one(T)
    ) where T
    α = convert(T, α)
    @inbounds for i in 1:min(size(m)...)
        m[i, i] += α
    end
    nothing
end


# =========================================================================================
"""
    isdegenerate(elem::Triangle{T})

Tests whether the given triangle is degenerate.

# Return type
`Bool`
"""
function isdegenerate(elem::Triangle{T}) where T
    @assert length(elem.v1) == length(elem.v2) == length(elem.v3) == 3
    u1 = elem.v2 .- elem.v1
    u2 = elem.v3 .- elem.v1
    cosine = _dot(u1, u2) / _norm(u1) / _norm(u2)
    _norm(elem.v1 .- elem.v2) < _etol(T) || _norm(elem.v1 .- elem.v3) < _etol(T) ||
    _norm(elem.v2 .- elem.v3) < _etol(T) || 1 - abs(cosine) <= _etol(T)
end


# =========================================================================================
"""
    _seek(
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
function _seek(fh::IOStream, prefix::String, skiptheline::Bool=true)
    m = -1
    found = false
    while !eof(fh)
        skiptheline || (m = position(fh))
        startswith(readline(fh), prefix) && (found = true) && break
    end
    found && (skiptheline || seek(fh, m))
    nothing
end


# =========================================================================================
"""
    _cos(
        u    ::AbstractVector{T},
        v    ::AbstractVector{T},
        unorm::T=_norm(u),
        vnorm::T=_norm(v)
    )

Computes the cosine of the angle between the given vectors `u` and `v` with lengths `unorm`
and `vnorm`, respectively.

# Return type
`T`
"""
@inline function _cos(
        u    ::AbstractVector{T},
        v    ::AbstractVector{T},
        unorm::T=_norm(u),
        vnorm::T=_norm(v)
    ) where T
    _dot(u, v) / (unorm * vnorm)
end


# =========================================================================================
"""
    cathetus(hyp::T, cosθ::T)

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
@inline function cathetus(hyp::T, cosθ::T) where T
    √(hyp^2 * (1 - cosθ^2))
end


# =========================================================================================
"""
    _sign(
        u::AbstractVector{T},
        v::AbstractVector{T},
        n::AbstractVector{T}
    )

Determines whether the normal vector of the plane specified by the vectors `u` and `v` has
the same orientation as the given normal vector `n`. Returns ``1`` if both normals have the
same orientation, ``0`` if at least one of the vectors is zero, and ``-1`` otherwise.

# Return type
`T`
"""
@inline function _sign(
        u::AbstractVector{T},
        v::AbstractVector{T},
        n::AbstractVector{T}
    ) where T
    # Devectorized version of sign((u1 × u2) ⋅ normal)
    sign(
        (u[2]*v[3] - u[3]*v[2]) * n[1] +
        (u[3]*v[1] - u[1]*v[3]) * n[2] +
        (u[1]*v[2] - u[2]*v[1]) * n[3]
    )
end


# =========================================================================================
"""
    distance(
        q   ::AbstractVector{T},
        elem::Triangle{T}
    )

Calculates the (positive or negative) distance from the given point `q` to the plane the
given triangle `elem` is located in. Negative distances correspond to point locations on
the 'inside' wrt. to the element's (outward-facing) normal vector.

# Return type
`T`
"""
@inline function distance(q::AbstractVector{T}, elem::Triangle{T}) where T
    _dot(q, elem.normal) - elem.distorig
end


# =========================================================================================
"""
    ddot(
        u::AbstractVector{T},
        v::AbstractVector{T},
        n::AbstractVector{T}
    )

Devectorized computation of `(u-v)⋅n`.

# Return type
`T`
"""
@inline function ddot(
        u::AbstractVector{T},
        v::AbstractVector{T},
        n::AbstractVector{T}
    ) where T
    (u[1] - v[1]) * n[1] + (u[2] - v[2]) * n[2] + (u[3] - v[3]) * n[3]
end


# =========================================================================================
"""
    _dot(
        u::AbstractVector{T},
        v::AbstractVector{T}
    )

Fast dot product for 3-element vectors.

# Return type
`T`
"""
@inline function _dot(u::AbstractVector{T}, v::AbstractVector{T}) where T
    u[1] * v[1] + u[2] * v[2] + u[3] * v[3]
end


# =========================================================================================
"""
    _norm(u::AbstractVector{T})

Fast Euclidean norm for 3-element vectors.

# Return type
`T`
"""
@inline function _norm(u::AbstractVector)
    √_dot(u, u)
end


# =========================================================================================
"""
    _reverseindex(v::AbstractVector{T})

Creates a reverse index for the given vector `v`, that is, a dictionary linking the object
IDs of the vector elements to the corresponding position in the vector.

# Return type
`IdDict`
"""
@inline function _reverseindex(v::AbstractVector)
    IdDict(n => i for (i, n) in enumerate(v))
end


# =========================================================================================
"""
    obspoints_line(
        u::AbstractVector{T},
        v::AbstractVector{T},
        n::Int
    )

Generates `n` evenly distributed observation points along the line segment from `u` to `v`.

# Return type
`Generator -> Vector{T}`

# Example
```julia
for ξ in obspoints_line([0, 0, 0], [1, 1, 1], 10)
    ...
end
```
"""
@inline function obspoints_line(u::AbstractVector{T}, v::AbstractVector{T}, n::Int) where T
    (u .+ T(i) .* (v .- u) for i in LinRange(0, 1, n))
end


# =========================================================================================
"""
    obspoints_plane(
        a  ::AbstractVector{T},
        b  ::AbstractVector{T},
        c  ::AbstractVector{T},
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
`Generator -> Vector{T}`

# Example
```julia
for Ξ in obspoints_plane(...)
    for ξ in Ξ
        ...
    end
end
```
"""
@inline function obspoints_plane(
        a  ::AbstractVector{T},
        b  ::AbstractVector{T},
        c  ::AbstractVector{T},
        nba::Int,
        nbc::Int
    ) where T
    (obspoints_line(ξ, c .+ (ξ .- b), nbc) for ξ in obspoints_line(b, a, nba))
end


# =========================================================================================
# Convenience aliases

const _axpy! = BLAS.axpy!

@inline _gemm(
    α::T,
    a::AbstractMatrix{T},
    b::AbstractMatrix{T},
) where T = BLAS.gemm('N', 'N', α, a, b)

@inline _gemv(
    α::T,
    m::AbstractMatrix{T},
    v::AbstractVector{T}
) where T = BLAS.gemv('N', α, m, v)

@inline _gemv!(
    α   ::T,
    m   ::AbstractMatrix{T},
    v   ::AbstractVector{T},
    β   ::T,
    dest::AbstractVector{T}
) where T = BLAS.gemv!('N', α, m, v, β, dest)

@inline _gemv!(
    α   ::T,
    m   ::AbstractMatrix{T},
    v   ::AbstractVector{T},
    dest::AbstractVector{T}
) where T = _gemv!(α, m, v, one(T), dest)

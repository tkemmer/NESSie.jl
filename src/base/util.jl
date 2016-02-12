#=
    Computes the element properties, that is, centroid, normal, distance to origin, and
    area.

    @param elem
                Triangle of interest
=#
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

#=
    Returns a dictionary that links each node pointer to the corresponding index in
    the `nodes` list. Registers only the last index if the same node is included
    multiple times.

    @param nodes
        List of nodes
    @return Dict{Pointer{T}, Int}
=#
indexmap{T}(nodes::Vector{Vector{T}}) = Dict([pointer(node) => i for (i, node) in enumerate(nodes)])

#=
    Unpacks the given vector of vectors into a single vector:
    [[1, 2, 3], [4, 5, 6]] => [1, 2, 3, 4, 5, 6]

    @param data
        Vector of vectors
    @param innerdim
        Dimension of the inner vectors
    @return Vector{T}
=#
unpack{T}(data::Vector{Vector{T}}, innerdim=3) = T[o for o in T[o[i] for i in 1:innerdim, o in data]]

#=
    Returns a list containing the normal vectors of the given nodes with respect to the
    given surface elements.

    @param nodes
        List of nodes
    @param elements
        List of surface elements
    @param invert
        Specifies whether all normals should be inverted
    @return
        Vector{Vector{T}}
=#
function vertexnormals{T}(nodes::Vector{Vector{T}}, elements::Vector{Triangle{T}}, invert::Bool=false)
    revidx = reverseindex(nodes)
    normals = Vector{T}[zeros(T, 3) for _ in 1:length(nodes)]
    count = zeros(T, length(nodes))
    @inbounds for elem in elements, node in (elem.v1, elem.v2, elem.v3)
        idx = revidx[object_id(node)]
        count[idx] += 1
        for i in 1:3
            normals[idx][i] += (elem.normal[i]-normals[idx][i]) / count[idx]
        end
    end
    invert ? -normals : normals
end

#=
    Initializes the given matrix m with αI, with I being an identity
    matrix with the same dimensions as m.

    @param m
        Corresponding matrix
    @param α
        Coefficient of the identity matrix
=#
function eye!{T}(m::Union{DenseArray{T,2}, SubArray{T,2}}, α::Number=one(T))
    fill!(m, zero(T))
    α = convert(T, α)
    @inbounds for i in 1:min(size(m)...)
        m[i, i] = α
    end
    nothing
end

#=
    Tests whether the triangle with the given nodes is degenerate.

    @param v1
        First node of the triangle
    @param v2
        Second node of the triangle
    @param v3
        Third node of the triangle
    @return bool
=#
function isdegenerate{T <: AbstractFloat}(v1::Vector{T}, v2::Vector{T}, v3::Vector{T})
    @assert length(v1) == length(v2) == length(v3) == 3
    u1 = v2 - v1
    u2 = v3 - v1
    cosine = u1 ⋅ u2 / vecnorm(u1) / vecnorm(u2)
    v1 == v2 || v1 == v3 || v2 == v3 || 1 - abs(cosine) <= 1e-10
end
isdegenerate{T}(elem::Triangle{T}) = isdegenerate(elem.v1, elem.v2, elem.v3)

#=
    Fast-forwards an IOStream to the next line starting with the
    given prefix. In case there is no such line. the stream handle
    will be set to EOF.

    @param fh
        A stream object
    @param
        Prefix of interest
    @param skiptheline
        If true, said line will also be skipped
=#
function seek(fh::IOStream, prefix::ASCIIString, skiptheline::Bool=true)
    m = -1
    while !eof(fh)
        skiptheline || (m = position(fh))
        startswith(readline(fh), prefix) && break
    end
    eof(fh) || skiptheline || seek(fh, m)
    nothing
end

#=
    Computes the cosine of the angle between the given vectors.

    @param u
        First vector
    @param unorm
        Length of the first vector
    @param v
        Second vector
    @param vnorm
        Length of the second vector
    @return T
=#
cos{T}(u::Vector{T}, unorm::T, v::Vector{T}, vnorm::T) = (u ⋅ v / (unorm * vnorm))
cos{T}(u::Vector{T}, v::Vector{T}) = cos(u, vecnorm(u), v, vecnorm(v))

#=
    Computes the cathetus c1 of a triangle given the hypotenuse h and the cosine of the
    exterior angle θ between the hypotenuse and the other cathetus c2.

    c2 = h * cos θ
    h² = c1² + c2²
    <=> c1 = √(h² * (1 - cos² θ))

    @param hyp
        Hypotenuse of the triangle
    @param cosθ
        Cosine of the exterior angle between the hypotenuse and the other cathetus
    @return T
=#
cathetus{T}(hyp::T, cosθ::T) = √(hyp^2 * (1 - cosθ^2))

#=
    Determines whether the normal vector of the plane specified by the vectors u1 and u2
    has the same orientation as the given normal vector.

    @param u1
        First vector
    @param u2
        Second vector
    @param
        Normal vector to compare against
    @return 1 if both normals have the same orientation, 0 if at least one of the
        vectors is zero, -1 otherwise.
=#
sign{T}(u1::Vector{T}, u2::Vector{T}, normal::Vector{T}) = sign((u1 × u2) ⋅ normal)

#=
    Calculates the (positive or negative) distance from the given point q to the plane
    given in Hesse normal form, that is, in the form of a unit normal vector and its
    distance to the origin.

    @param q
        Point of interest
    @param normal
        Unit normal vector of the plane
    @param distorig
        Distance from the origin to the plane (≥ 0)
    @return T
=#
distance{T}(q::Vector{T}, normal::Vector{T}, distorig::T) = q ⋅ normal - distorig
distance{T}(q::Vector{T}, elem::Triangle{T}) = distance(q, elem.normal, elem.distorig)

#=
    Creates a reverse index for the given vector, that is, a dictionary linking the
    object IDs of the vector elements to the corresponding position in the vector.

    @param v
        A vector
    @return Dict{UInt, UInt}
=#
reverseindex{T}(v::Vector{T}) = Dict{UInt, UInt}([object_id(e) => i for (i,e) in enumerate(v)])

# Convenience aliases
gemv!{T}(α::T, m::Union{DenseArray{T,2}, SubArray{T,2}}, v::Vector{T}, dest::Union{DenseArray{T,1}, SubArray{T,1}}) = gemv!(α, m, v, one(T), dest)
gemv!{T}(α::T, m::Union{DenseArray{T,2}, SubArray{T,2}}, v::Vector{T}, β::T, dest::Union{DenseArray{T,1}, SubArray{T,1}}) = gemv!('N', α, m, v, β, dest)
# =========================================================================================
"""
    readstl(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads a surface model from the given STL file.

!!! note
    This file type does not support charge models! Hence, the charge list of the returning `Model` object is empty and has to be set separately.

# Specification
<https://fabbers.com/tech/STL_Format>

# Return type
`Model{T, Triangle{T}}`

# Alias

    readstl(
        fname::String,
             ::Type{T}=Float64
    )

Reads the model using a file name rather than an `IOStream` object.
"""
function readstl(
    stream::IOStream,
          ::Type{T}=Float64
    ) where T
    # skip header
    skip(stream, 80)
    # read number of triangles
    numelem = read(stream, UInt32)
    elements = Vector{Triangle{T}}(undef, numelem)
    nodes = Set{Vector{T}}()

    for i in 1:numelem
        skip(stream, 12) # skip normal (will be computed by props)
        v1 = map(T, [read(stream, Float32) for _ in 1:3])
        v2 = map(T, [read(stream, Float32) for _ in 1:3])
        v3 = map(T, [read(stream, Float32) for _ in 1:3])
        skip(stream, 2)

        # make sure nodes with the same coordinates refer to the same object
        v1 ∈ nodes || push!(nodes, v1)
        v2 ∈ nodes || push!(nodes, v2)
        v3 ∈ nodes || push!(nodes, v3)
        elements[i] = Triangle(
            getkey(nodes.dict, v1, v1),
            getkey(nodes.dict, v2, v2),
            getkey(nodes.dict, v3, v3)
        )
    end

    Model(collect(Vector{T}, nodes), elements)
end

@inline function readstl(
    fname::String,
         ::Type{T}=Float64
) where T
    open(fh -> readstl(fh, T), fname)
end


# =========================================================================================
"""
    writestl(
        stream::IOStream,
        model ::Model{T, Triangle{T}}
    )

Creates a binary STL file from a given surface model.

# Specification
<http://www.fabbers.com/tech/STL_Format>

# Return type
`Void`

# Alias

    writestl(
        fname::String,
        model::Model{T, Triangle{T}}
    )

Creates the STL file by name rather than `IOStream` object.
"""
function writestl(
    stream::IOStream,
    model ::Model{T, Triangle{T}}
) where T
    # header (empty)
    write(stream, zeros(UInt8, 80))
    # number of triangles
    write(stream, UInt32(length(model.elements)))

    for elem in model.elements
        write(stream, map(Float32, elem.normal))
        write(stream, map(Float32, elem.v1))
        write(stream, map(Float32, elem.v2))
        write(stream, map(Float32, elem.v3))
        write(stream, zero(UInt16))
    end
    nothing
end

@inline function writestl(fname::String, model::Model{T, Triangle{T}}) where T
    open(fh -> writestl(fh, model), fname, "w")
end

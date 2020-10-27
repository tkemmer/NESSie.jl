# =========================================================================================
"""
    readstl{T}(
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

    readstl{T}(
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

    for i in 1:numelem
        skip(stream, 12) # skip normal (will be computed by props)
        elements[i] = Triangle(
            map(T, [read(stream, Float32) for _ in 1:3]),
            map(T, [read(stream, Float32) for _ in 1:3]),
            map(T, [read(stream, Float32) for _ in 1:3])
        )
        skip(stream, 2)
    end

    nodes = collect(Set(unpack([[e.v1, e.v2, e.v3] for e in elements])))
    Model(nodes, elements)
end

function readstl(
    fname::String,
         ::Type{T}=Float64
) where T
    open(fh -> readstl(fh, T), fname)
end


# =========================================================================================
"""
    writestl{T}(
        stream::IOStream,
        model ::Model{T, Triangle{T}}
    )

Creates a binary STL file from a given surface model.

# Specification
<http://www.fabbers.com/tech/STL_Format>

# Return type
`Void`

# Alias

    writestl{T}(
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

function writestl(fname::String, model::Model{T, Triangle{T}}) where T
    open(fh -> writestl(fh, model), fname, "w")
end

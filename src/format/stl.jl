# =========================================================================================
"""
    function writestl{T}(
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

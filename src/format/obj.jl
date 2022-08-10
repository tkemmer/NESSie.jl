# =========================================================================================
"""
    function writeobj{T}(
        stream::IOStream,
        model ::Model{T, Triangle{T}}
    )

Creates a Wavefront OBJ file from a given surface model.

# Specification
    <https://www.loc.gov/preservation/digital/formats/fdd/fdd000507.shtml>
    <http://paulbourke.net/dataformats/obj/>

# Return type
`Void`

# Alias

    writeobj{T}(
        fname::String,
        model::Model{T, Triangle{T}}
    )

Creates the OBJ file by name rather than `IOStream` object.
"""
function writeobj(
    stream::IOStream,
    model ::Model{T, Triangle{T}}
) where T
    # nodes
    for node in model.nodes
        println(stream, "v ", node[1], " ", node[2], " ", node[3])
    end
    # vertex normals
    vns = vertexnormals(model)
    for vn in vns
        println(stream, "vn ", vn[1], " ", vn[2], " ", vn[3])
    end
    # elements
    ridx = reverseindex(model.nodes)
    for elm in model.elements
        print(stream, "f ")
        print(stream, ridx[objectid(elm.v1)], "//", ridx[objectid(elm.v1)], " ")
        print(stream, ridx[objectid(elm.v2)], "//", ridx[objectid(elm.v2)], " ")
        println(stream, ridx[objectid(elm.v3)], "//", ridx[objectid(elm.v3)])
    end
end

@inline function writeobj(fname::String, model::Model{T, Triangle{T}}) where T
    open(fh -> writeobj(fh, model), fname, "w")
end

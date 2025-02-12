# =========================================================================================
"""
    writeskel(
        stream::IOStream,
        model ::Model{T}
    )

Creates a SKEL file from a given surface or volume model, representing the model as a
collection of points and polylines.

# Specification
<http://www.geomview.org/docs/html/SKEL.html>

# Return type
`Nothing`

# Alias

    writeskel(
        fname::String,
        model::Model{T}
    )

Creates the SKEL file by name rather than `IOStream` object.
"""
function writeskel(
        stream::IOStream,
        model ::Model{T, Triangle{T}}
    ) where T
    println(stream, "SKEL")
    println(stream, "$(length(model.nodes))\t$(length(model.elements))")
    for node in model.nodes
        println(stream, "$(node[1])\t$(node[2])\t$(node[3])")
    end
    revidx = reverseindex(model.nodes)
    for elem in model.elements
        (v1, v2, v3) = (revidx[objectid(elem.v1)] - 1,
                        revidx[objectid(elem.v2)] - 1,
                        revidx[objectid(elem.v3)] - 1)
        println(stream, "4\t$v1\t$v2\t$v3\t$v1")
    end
    nothing
end

function writeskel(
        stream::IOStream,
        model ::Model{T, Tetrahedron{T}}
    ) where T
    println(stream, "SKEL")
    println(stream, "$(length(model.nodes))\t$(length(model.elements))")
    for node in model.nodes
        println(stream, "$(node[1])\t$(node[2])\t$(node[3])")
    end
    revidx = reverseindex(model.nodes)
    for elem in model.elements
        (v1, v2, v3, v4) = (revidx[objectid(elem.v1)] - 1,
                            revidx[objectid(elem.v2)] - 1,
                            revidx[objectid(elem.v3)] - 1,
                            revidx[objectid(elem.v4)] - 1)
        println(stream, "7\t$v1\t$v2\t$v3\t$v4\t$v2\t$v4\t$v1")
    end
    nothing
end

@inline function writeskel(
        fname::String,
        model::M
    ) where {T, M <: Model{T}}
    open(fh -> writeskel(fh, model), fname, "w")
end

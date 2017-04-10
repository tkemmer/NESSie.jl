#=
    Creates a SKEL file from a given surface model, representing the model as a collection of points and polylines.

    File specification:
    http://www.geomview.org/docs/html/SKEL.html

    @param fname/stream
        Path or handle to (writable) SKEL file
    @param model
        A surface model
    @return nothing
=#
function writeskel{T}(stream::IOStream, model::SurfaceModel{T})
    println(stream, "SKEL")
    println(stream, "$(length(model.nodes))\t$(length(model.elements))")
    for node in model.nodes
        println(stream, "$(node[1])\t$(node[2])\t$(node[3])")
    end
    revidx = reverseindex(model.nodes)
    for elem in model.elements
        (v1, v2, v3) = (revidx[object_id(elem.v1)] - 1,
                        revidx[object_id(elem.v2)] - 1,
                        revidx[object_id(elem.v3)] - 1)
        println(stream, "4\t$v1\t$v2\t$v3\t$v1")
    end
    nothing
end
writeskel{T}(fname::String, model::SurfaceModel{T}) = open(fh -> writeskel(fh, model), fname, "w")

#=
    Creates a SKEL file from a given volume model, representing the model as a collection of points and polylines.

    File specification:
    http://www.geomview.org/docs/html/SKEL.html

    @param fname/stream
        Path or handle to (writable) SKEL file
    @param model
        A volume model
    @return nothing
=#
function writeskel{T}(stream::IOStream, model::VolumeModel{T})
    println(stream, "SKEL")
    println(stream, "$(length(model.nodes))\t$(length(model.elements))")
    for node in model.nodes
        println(stream, "$(node[1])\t$(node[2])\t$(node[3])")
    end
    revidx = reverseindex(model.nodes)
    for elem in model.elements
        (v1, v2, v3, v4) = (revidx[object_id(elem.v1)] - 1,
                            revidx[object_id(elem.v2)] - 1,
                            revidx[object_id(elem.v3)] - 1,
                            revidx[object_id(elem.v4)] - 1)
        println(stream, "7\t$v1\t$v2\t$v3\t$v4\t$v2\t$v4\t$v1")
    end
    nothing
end
writeskel{T}(fname::String, model::VolumeModel{T}) = open(fh -> writeskel(fh, model), fname, "w")

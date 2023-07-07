# =========================================================================================
"""
    readoff{T <: AbstractFloat}(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads a surface model from the given OFF file.

!!! note
    This file type does not support charge models! Hence, the charge list of the returning
    `Model` object is empty and has to be set separately.

# Specification
<http://www.geomview.org/docs/html/OFF.html>

# Return type
`Model{T}`

# Alias

    readoff{T}(fname::String, ::Type{T}=Float64)

Reads the model using a file name rather than a `IOStream` object.
"""
function readoff(
        stream::IOStream,
              ::Type{T}=Float64
    ) where T <: AbstractFloat
    eof(stream) && return Model{T, Triangle{T}}()
    @assert readline(stream) == "OFF" "Invalid OFF file"

    # read number of nodes and elements
    numnodes, numelem = [parse(Int, s) for s in split(readline(stream))]

    nodes = readoff_nodes(stream, numnodes, T)
    Model(nodes, readoff_elements(stream, numelem, nodes))
end

function readoff(
        fname::String,
             ::Type{T}=Float64
    ) where T
    open(fh -> readoff(fh, T), fname)
end


# =========================================================================================
"""
    readoff_nodes{T <: AbstractFloat}(
        stream::IOStream,
        n     ::Int,
              ::Type{T}=Float64
    )

Reads the first `n` nodes from the given OFF file.

# Return type
`Vector{Vector{T}}`
"""
function readoff_nodes(
        stream::IOStream,
        n     ::Int,
              ::Type{T}=Float64
    ) where T <: AbstractFloat
    Vector{T}[[parse(T, e) for e in split(readline(stream))] for _ in 1:n]
end


# =========================================================================================
"""
    readoff_elements{T <: AbstractFloat}(
        stream::IOStream,
        n     ::Int,
        nodes ::Vector{Vector{T}},
              ::Type{T}=Float64
    )

Reads the first `n` elements from the given OFF file.

# Return type
`Vector{Triangle{T}}`
"""
function readoff_elements(
        stream::IOStream,
        n     ::Int,
        nodes ::Vector{Vector{T}}
    ) where T <: AbstractFloat
    Triangle{T}[
        Triangle(
            [nodes[parse(Int, e) + 1] for e in split(readline(stream))[2:end]]...
        ) for _ in 1:n
    ]
end

#=
    TODO
=#
function writeoff{T}(stream::IOStream, model::SurfaceModel{T})

end

#=
    TODO
=#
function writeoff{T}(stream::IOStream, model::VolumeModel{T})
    println(stream, "OFF")
    println(stream, "$(length(model.nodes))\t$(length(model.elements) * 4)\t$(length(model.elements) * 12)")
    for node in model.nodes
        println(stream, "$(node[1])\t$(node[2])\t$(node[3])")
    end
    revidx = reverseindex(model.nodes)
    for elem in model.elements
        (v1, v2, v3, v4) = (revidx[object_id(elem.v1)] - 1,
                            revidx[object_id(elem.v2)] - 1,
                            revidx[object_id(elem.v3)] - 1,
                            revidx[object_id(elem.v4)] - 1)
        println(stream, "3\t$v1\t$v3\t$v2")
        println(stream, "3\t$v1\t$v2\t$v4")
        println(stream, "3\t$v2\t$v3\t$v4")
        println(stream, "3\t$v3\t$v1\t$v4")
    end
    nothing
end
writeoff{T}(fname::String, model::VolumeModel{T}) = open(fh -> writeoff(fh, model), fname, "w")


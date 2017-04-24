# =========================================================================================
"""
    readmsms{T <: AbstractFloat}(
        vertstream::IOStream,
        facestream::IOStream,
                  ::Type{T}=Float64
    )

Reads a surface model from the given MSMS-generated `.face` and `.vert` files.

!!! note
    This file type does not support charge models! Hence, the charge list of the returning
    `SurfaceModel` object is empty and has to be set separately.

# Return type
`SurfaceModel{T}`

# Alias

    readmsms{T}(fname::String, ::Type{T}=Float64)

Reads the model using a common file name prefix (`fname.{vert,face}`) for both files
rather than `IOStream` objects.

"""
function readmsms{T <: AbstractFloat}(
        vertstream::IOStream,
        facestream::IOStream,
        ::Type{T}=Float64
    )
    nodes = readmsms_nodes(vertstream, T)
    SurfaceModel(nodes, readmsms_elements(facestream, nodes, T), Charge{T}[])
end

function readmsms{T <: AbstractFloat}(fname::String, ::Type{T}=Float64)
    open(ff -> open(fv -> readmsms(fv, ff, T), "$fname.vert"), "$fname.face")
end


# =========================================================================================
"""
    readmsms_nodes{T <: AbstractFloat}(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads all nodes from the given MSMS-generated `.vert` file.

# Return type
`Vector{Vector{T}}`
"""
function readmsms_nodes{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    nodes = Vector{T}[]
    # skip header lines
    [readline(stream) for i in 1:3]
    for line in eachline(stream)
        push!(nodes, [parse(T, e) for e in split(line)[1:3]])
    end
    nodes
end


# =========================================================================================
"""
    readmsms_elements{T <: AbstractFloat}(
        stream::IOStream,
        nodes ::Vector{Vector{T}},
              ::Type{T}=Float64
    )

Reads all elements from the given MSMS-generated `.face` file.

# Return type
`Vector{Triangle{T}}`
"""
function readmsms_elements{T <: AbstractFloat}(
        stream::IOStream,
        nodes::Vector{Vector{T}},
        ::Type{T}=Float64
    )
    elements = Triangle{T}[]
    # skip header lines
    [readline(stream) for i in 1:3]
    for line in eachline(stream)
        push!(elements, Triangle([nodes[parse(Int, e)] for e in split(line)[1:3]]...))
    end
    elements
end

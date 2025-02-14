# =========================================================================================
"""
    readmsms(
        vertstream::IOStream,
        facestream::IOStream,
                  ::Type{T}=Float64
    )

Reads a surface model from the given MSMS-generated `.face` and `.vert` files.

!!! note
    This file type does not support charge models! Hence, the charge list of the returning
    `Model` object is empty and has to be set separately.

# Return type
[`Model{T, Triangle{T}}`](@ref)

# Alias

    readmsms(fname::String, ::Type{T}=Float64)

Reads the model using a common file name prefix (`fname.{vert,face}`) for both files
rather than `IOStream` objects.

"""
function readmsms(
        vertstream::IOStream,
        facestream::IOStream,
                  ::Type{T}=Float64
    ) where T <: AbstractFloat
    nodes = readmsms_nodes(vertstream, T)
    Model(nodes, readmsms_elements(facestream, nodes))
end

@inline function readmsms(
        fname::String,
             ::Type{T}=Float64
    ) where T <: AbstractFloat
    open(ff -> open(fv -> readmsms(fv, ff, T), "$fname.vert"), "$fname.face")
end


# =========================================================================================
"""
    readmsms_nodes(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads all nodes from the given MSMS-generated `.vert` file.

# Return type
`Vector{Vector{T}}`
"""
function readmsms_nodes(
        stream::IOStream,
              ::Type{T}=Float64
    ) where T <: AbstractFloat
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
    readmsms_elements(
        stream::IOStream,
        nodes ::Vector{Vector{T}}
    )

Reads all elements from the given MSMS-generated `.face` file.

# Return type
[`Vector{Triangle{T}}`](@ref Triangle)
"""
function readmsms_elements(
        stream::IOStream,
        nodes ::Vector{Vector{T}}
    ) where T <: AbstractFloat
    elements = Triangle{T}[]
    # skip header lines
    [readline(stream) for i in 1:3]
    for line in eachline(stream)
        push!(elements, Triangle([nodes[parse(Int, e)] for e in split(line)[1:3]]...))
    end
    elements
end

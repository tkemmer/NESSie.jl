# =========================================================================================
"""
    readoff{T <: AbstractFloat}(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads a surface model from the given OFF file.

!!! note
    This file type does not support charge models! Hence, the charge list of the returning
    `SurfaceModel` object is empty and has to be set separately.

# Specification
<http://www.geomview.org/docs/html/OFF.html>

# Return type
`SurfaceModel{T}`

# Alias

    readoff{T}(fname::String, ::Type{T}=Float64)

Reads the model using a file name rather than a `IOStream` object.
"""
function readoff{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    eof(stream) && return SurfaceModel(Vector{T}[], Triangle{T}[], Charge{T}[])
    @assert readline(stream) == "OFF" "Invalid OFF file"

    # read number of nodes and elements
    numnodes, numelem = [parse(Int, s) for s in split(readline(stream))]

    nodes = readoff_nodes(stream, numnodes, T)
    SurfaceModel(nodes, readoff_elements(stream, numelem, nodes, T), Charge{T}[])
end
readoff{T}(fname::String, ::Type{T}=Float64) = open(fh -> readoff(fh, T), fname)


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
function readoff_nodes{T <: AbstractFloat}(stream::IOStream, n::Int, ::Type{T}=Float64)
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
function readoff_elements{T <: AbstractFloat}(
        stream::IOStream,
        n::Int,
        nodes::Vector{Vector{T}},
        ::Type{T}=Float64
    )
    Triangle{T}[
        Triangle(
            [nodes[parse(Int, e) + 1] for e in split(readline(stream))[2:end]]...
        ) for _ in 1:n
    ]
end

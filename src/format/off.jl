# =========================================================================================
"""
    readoff(
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
[`Model{T}`](@ref)

# Alias

    readoff(fname::String, ::Type{T}=Float64)

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

@inline function readoff(
        fname::String,
             ::Type{T}=Float64
    ) where T
    open(fh -> readoff(fh, T), fname)
end


# =========================================================================================
"""
    readoff_nodes(
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
    readoff_elements(
        stream::IOStream,
        n     ::Int,
        nodes ::Vector{Vector{T}},
              ::Type{T}=Float64
    )

Reads the first `n` elements from the given OFF file.

# Return type
[`Vector{Triangle{T}}`](@ref Triangle)
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

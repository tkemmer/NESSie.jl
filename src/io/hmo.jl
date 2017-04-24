# =========================================================================================
"""
    readhmo{T <: AbstractFloat}(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads a complete surface model from the given HMO file.

# Return type
`SurfaceModel{T}`

# Alias

    readhmo{T}(fname::String, ::Type{T}=Float64)

Reads the model using a file name rather than a `IOStream` object.
"""
function readhmo{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    nodes = readhmo_nodes(stream, T)
    SurfaceModel(nodes, readhmo_elements(stream, nodes, T), readhmo_charges(stream, T))
end
readhmo{T}(fname::String, ::Type{T}=Float64) = open(fh -> readhmo(fh, T), fname)


# =========================================================================================
"""
    readhmo_nodes{T <: AbstractFloat}(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads all nodes from the given HMO file.

# Return type
`Vector{Vector{T}}`
"""
function readhmo_nodes{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    nodes = Vector{T}[]
    seek(stream, "BEG_NODL_DATA")
    readline(stream) # skip first line
    for line in eachline(stream)
        line == "END_NODL_DATA" && break
        push!(nodes, [parse(T, a) for a in split(line)[2:end]])
    end
    nodes
end


# =========================================================================================
"""
    readhmo_elements{T <: AbstractFloat}(
        stream::IOStream,
        nodes ::Vector{Vector{T}},
              ::Type{T}=Float64
    )

Reads all elements from the given HMO file.

# Return type
`Vector{Triangle{T}}`
"""
function readhmo_elements{T <: AbstractFloat}(
        stream::IOStream,
        nodes::Vector{Vector{T}},
        ::Type{T}=T
    )
    elements = Triangle{T}[]
    seek(stream, "BEG_ELEM_DATA")
    readline(stream) # skip first line
    for line in eachline(stream)
        line == "END_ELEM_DATA" && break
        push!(elements, Triangle([nodes[parse(Int, a)] for a in split(line)[4:end]]...))
    end
    elements
end


# =========================================================================================
"""
    readhmo_charges{T <: AbstractFloat}(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads all charges from the given HMO file.

# Return type
`Vector{Charge{T}}`
"""
function readhmo_charges{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    charges = Charge{T}[]
    seek(stream, "BEG_CHARGE_DATA")
    readline(stream) # skip first line
    for line in eachline(stream)
        line == "END_CHARGE_DATA" && break
        push!(charges, Charge([parse(T, a) for a in split(line)[2:end]]...))
    end
    charges
end

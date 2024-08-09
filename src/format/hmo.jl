# =========================================================================================
"""
    readhmo(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads a complete surface model from the given HMO file.

# Return type
`Model{T, Triangle{T}}`

# Alias

    readhmo(fname::String, ::Type{T}=Float64)

Reads the model using a file name rather than a `IOStream` object.
"""
function readhmo(
        stream::IOStream,
              ::Type{T}=Float64
    ) where T <: AbstractFloat
    nodes = readhmo_nodes(stream, T)
    Model(nodes, readhmo_elements(stream, nodes), readhmo_charges(stream, T))
end

@inline function readhmo(
        fname::String,
             ::Type{T}=Float64
    ) where T
    open(fh -> readhmo(fh, T), fname)
end


# =========================================================================================
"""
    readhmo_nodes(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads all nodes from the given HMO file.

# Return type
`Vector{Vector{T}}`
"""
function readhmo_nodes(
        stream::IOStream,
              ::Type{T}=Float64
    ) where T <: AbstractFloat
    nodes = Vector{T}[]
    _seek(stream, "BEG_NODL_DATA")
    readline(stream) # skip first line
    for line in eachline(stream)
        line == "END_NODL_DATA" && break
        push!(nodes, [parse(T, a) for a in split(line)[2:end]])
    end
    nodes
end


# =========================================================================================
"""
    readhmo_elements(
        stream::IOStream,
        nodes ::Vector{Vector{T}}
    )

Reads all elements from the given HMO file.

# Return type
`Vector{Triangle{T}}`
"""
function readhmo_elements(
        stream::IOStream,
        nodes ::Vector{Vector{T}}
    ) where T <: AbstractFloat
    elements = Triangle{T}[]
    _seek(stream, "BEG_ELEM_DATA")
    readline(stream) # skip first line
    for line in eachline(stream)
        line == "END_ELEM_DATA" && break
        push!(elements, Triangle([nodes[parse(Int, a)] for a in split(line)[4:end]]...))
    end
    elements
end


# =========================================================================================
"""
    readhmo_charges(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads all charges from the given HMO file.

# Return type
`Vector{Charge{T}}`
"""
function readhmo_charges(
        stream::IOStream,
              ::Type{T}=Float64
    ) where T <: AbstractFloat
    charges = Charge{T}[]
    _seek(stream, "BEG_CHARGE_DATA")
    readline(stream) # skip first line
    for line in eachline(stream)
        line == "END_CHARGE_DATA" && break
        push!(charges, Charge([parse(T, a) for a in split(line)[2:end]]...))
    end
    charges
end


# =========================================================================================
"""
    writehmo(
        stream::IOStream,
        model ::Model{T, Triangle{T}}
    )

Creates a HMO file from the given surface model.

# Return type
`Void`

# Alias

    writehmo(
        fname::String,
        model::Model{T, Triangle{T}}
    )

Creates the HMO file by name rather than `IOStream` object.
"""
function writehmo(
        stream::IOStream,
        model ::Model{T, Triangle{T}}
    ) where T
    write(stream, "BEG_NODL_DATA\n")
    write(stream, "$(length(model.nodes))\n")
    for (i, node) in enumerate(model.nodes)
        join(stream, Any[i, node...], "\t")
        write(stream, "\n")
    end
    write(stream, "END_NODL_DATA\n\n")

    getidx = v -> findfirst(isequal(v), model.nodes)

    write(stream, "BEG_ELEM_DATA\n")
    join(stream, [length(model.elements), zeros(Int, 12)...], "\t")
    write(stream, "\n")
    for (i, elem) in enumerate(model.elements)
        join(stream, Any[i, 1, 103, getidx(elem.v1), getidx(elem.v2), getidx(elem.v3)], "\t")
        write(stream, "\n")
    end
    write(stream, "END_ELEM_DATA\n\n")

    write(stream, "BEG_CHARGE_DATA\n")
    write(stream, "$(length(model.charges))\n")
    for (i, charge) in enumerate(model.charges)
        join(stream, Any[i, charge.pos..., charge.val], "\t")
        write(stream, "\n")
    end
    write(stream, "END_CHARGE_DATA\n")
    return
end

@inline function writehmo(fname::String, model::Model{T, Triangle{T}}) where T
    open(fh -> writehmo(fh, model), fname, "w")
end

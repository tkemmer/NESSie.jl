# =========================================================================================
"""
    readmcsf{T <: AbstractFloat}(
        stream::IOStream,
              ::Type{T}=Float64;
        # kwargs
        domain::Symbol=:none
    )

Reads a volume model from the given GAMer-generated mcsf file.

!!! note
    This file type does not support charge models! Hence, the charge list of the returning
    `VolumeModel` object is empty and has to be set separately.

# Arguments
 * `domain` Element domain (solute `:Ω`, solvent `:Σ`, or `:none`)

# Return type
`VolumeModel{T}`

# Aliases

    readmcsf{T}(
        fname ::String,
              ::Type{T}=Float64;
        # kwargs
        domain::Symbol=:none
    )

Reads the model using a file name rather than a `IOStream` object.

    readmcsf{T}(
        fnameΩ::String,
        fnameΣ::String,
        ::Type{T}=Float64
    )

Reads the contents of two separate files (given by name), sets the element domains to `:Ω`
or `:Σ`, respectively, and returns a single (merged) `VolumeModel` object.
"""
function readmcsf{T <: AbstractFloat}(
        stream::IOStream,
        ::Type{T}=Float64;
        domain::Symbol=:none
    )
    nodes = readmcsf_nodes(stream, T)
    VolumeModel(nodes, readmcsf_elements(stream, nodes, T, domain=domain), Charge{T}[])
end

function readmcsf{T}(
        fname::String,
        ::Type{T}=Float64;
        domain::Symbol=:none
    )
    open(fh -> readmcsf(fh, T, domain=domain), fname)
end

function readmcsf{T}(fnameΩ::String, fnameΣ::String, ::Type{T}=Float64)
    meshunion(readmcsf(fnameΩ, T, domain=:Ω), readmcsf(fnameΣ, T, domain=:Σ))
end


# =========================================================================================
"""
    readmcsf_nodes{T <: AbstractFloat}(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads all nodes from the given GAMer-generated mcsf file.

# Return type
`Vector{Vector{T}}`
"""
function readmcsf_nodes{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    nodes = Vector{T}[]
    seek(stream, "vert=[")
    for line in eachline(stream)
        startswith(line, "%") && continue        # skip comments
        startswith(line, "];") && break          # all nodes read
        push!(nodes, [parse(T, e) for e in split(line)[3:5]])
    end
    nodes
end


# =========================================================================================
"""
    readmcsf_elements{T <: AbstractFloat}(
        stream::IOStream,
        nodes ::Vector{Vector{T}},
              ::Type{T}=Float64;
        # kwargs
        domain::Symbol=:none
    )

Reads all elements from the given GAMer-generated mcsf file.

# Return type
`Vector{Tetrahedron{T}}`
"""
function readmcsf_elements{T <: AbstractFloat}(
        stream::IOStream,
        nodes::Vector{Vector{T}},
        ::Type{T}=Float64;
        domain::Symbol=:none
    )
    elements = Tetrahedron{T}[]
    seek(stream, "simp=[")
    for line in eachline(stream)
        startswith(line, "%") && continue        # skip comments
        startswith(line, "];") && break          # all simplices read
        push!(
            elements,
            Tetrahedron([nodes[parse(Int, e)+1] for e in split(line)[8:end]]..., domain)
        )
    end
    elements
end

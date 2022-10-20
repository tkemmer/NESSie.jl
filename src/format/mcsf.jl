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
    `Model` object is empty and has to be set separately.

# Arguments
 * `domain` Element domain (solute `:Ω`, solvent `:Σ`, or `:none`)

# Return type
`Model{T, Tetrahedron{T}}`

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
or `:Σ`, respectively, and returns a single (merged) `Model` object.
"""
function readmcsf(
        stream::IOStream,
              ::Type{T}=Float64;
        domain::Symbol=:none
    ) where T <: AbstractFloat
    nodes = readmcsf_nodes(stream, T)
    Model(nodes, readmcsf_elements(stream, nodes, domain=domain))
end

@inline function readmcsf(
        fname ::String,
              ::Type{T}=Float64;
        domain::Symbol=:none
    ) where T
    open(fh -> readmcsf(fh, T, domain=domain), fname)
end

@inline function readmcsf(
        fnameΩ::String,
        fnameΣ::String,
              ::Type{T}=Float64
    ) where T
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
function readmcsf_nodes(
        stream::IOStream,
              ::Type{T}=Float64
    ) where T <: AbstractFloat
    nodes = Vector{T}[]
    _seek(stream, "vert=[")
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
        nodes ::Vector{Vector{T}};
        # kwargs
        domain::Symbol=:none
    )

Reads all elements from the given GAMer-generated mcsf file.

# Return type
`Vector{Tetrahedron{T}}`
"""
function readmcsf_elements(
        stream::IOStream,
        nodes ::Vector{Vector{T}};
        domain::Symbol=:none
    ) where T <: AbstractFloat
    elements = Tetrahedron{T}[]
    _seek(stream, "simp=[")
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

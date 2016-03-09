#=
    Reads the GAMer-generated volume mesh from the given mcsf file.

    @param stream
                Handle to mcsf file
    @param _
                Data type T for return value
    @return (Vector{Vector{T}}, Vector{Tetrahedron{T}})
=#
function readmcsf{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    nodes = readmcsf_nodes(stream, T)
    (nodes, readmcsf_elements(stream, nodes, T))
end

#=
    Reads the nodes of a GAMer-generated volume mesh from the given mcsf file.

    @param stream
                Handle to mcsf file
    @param _
                Data type T for return value
    @return Vector{Vector{T}}
=#
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

#=
    Reads the tetrahedra of a GAMer-generated volume mesh from the given mcsf file.

    @param stream
                Handle to mcsf file
    @param nodes
                List of reference nodes
    @param _
                Data type T for return value
    @return Vector{Tetrahedron{T}}
=#
function readmcsf_elements{T <: AbstractFloat}(stream::IOStream, nodes::Vector{Vector{T}}, ::Type{T}=Float64)
    elements = Tetrahedron{T}[]
    seek(stream, "simp=[")
    for line in eachline(stream)
        startswith(line, "%") && continue        # skip comments
        startswith(line, "];") && break          # all simplices read
        push!(elements, Tetrahedron([nodes[parse(Int, e)+1] for e in split(line)[8:end]]...))
    end
    elements
end

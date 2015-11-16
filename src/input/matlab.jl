#=
    Reads the GAMer-generated volume mesh from the given MATLAB file.

    @param stream
                Handle to MATLAB file
    @param _
                Data type T for return value
    @return (Vector{Vector{T}}, Vector{Tetrahedron{T}})
=#
function readmatlab{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    nodes = Vector{T}[]; elements = Tetrahedron{T}[]
    while !eof(stream)
        line = readline(stream)
        startswith(line, "vert=[") && (nodes = readmatlab_nodes(stream, T))
        startswith(line, "simp=[") && (elements = readmatlab_elements(stream, nodes, T))
    end
    (nodes, elements)
end

#=
    Reads the nodes of a GAMer-generated volume mesh from the given MATLAB file.

    @param stream
                Handle to MATLAB file
    @param _
                Data type T for return value
    @return Vector{Vector{T}}
=#
function readmatlab_nodes{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    nodes = Vector{T}[]
    while !eof(stream)
        line = readline(stream)
        startswith(line, "%") && continue
        startswith(line, "];") && break
        push!(nodes, [parse(T, e) for e in split(line)[3:5]])
    end
    nodes
end

#=
    Reads the tetrahedra of a GAMer-generated volume mesh from the given MATLAB file.

    @param stream
                Handle to MATLAB file
    @param nodes
                List of reference nodes
    @param _
                Data type T for return value
    @return Vector{Tetrahedron{T}}
=#
function readmatlab_elements{T <: AbstractFloat}(stream::IOStream, nodes::Vector{Vector{T}}, ::Type{T}=Float64)
    elements = Tetrahedron{T}[]
    while !eof(stream)
        line = readline(stream)
        startswith(line, "%") && continue
        startswith(line, "];") && break
        push!(elements, Tetrahedron([nodes[parse(Int, e)+1] for e in split(line)[8:end]]...))
    end
    elements
end

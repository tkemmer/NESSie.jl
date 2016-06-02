#=
    Reads the MSMS-generated surface mesh from the given .face and .vert files.

    @param vertstream
        Handle to MSMS .vert file
    @param facestream
        Handle to MSMS .face file
    @param _
        Data type T for return value
    @return (Vector{Vector{T}}, Vector{Tetrahedron{T}})
=#
function readmsms{T <: AbstractFloat}(vertstream::IOStream, facestream::IOStream, ::Type{T}=Float64)
    nodes = readmsms_nodes(vertstream, T)
    (nodes, readmsms_elements(facestream, nodes, T))
end

#=
    Convenience alias for readmsms(::IOStream, ::IOStream, ::Type{T}). Calls said
    method looking for the files ${fname}.vert and ${fname}.face.

    @param fname
        File name prefix of both .face and .vert files (${fname}.{vert,face})
    @param _
        Data type T for return value
=#
readmsms{T <: AbstractFloat}(fname::String, ::Type{T}=Float64) = open(ff -> open(fv -> readmsms(fv, ff, T), "$fname.vert"), "$fname.face")

#=
    Reads the nodes of a MSMS-generated surface mesh from the given .vert file.

    @param stream
        Handle to .vert file
    @param _
        Data type T for return value
    @return Vector{Vector{T}}
=#
function readmsms_nodes{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    nodes = Vector{T}[]
    # skip header lines
    [readline(stream) for i in 1:3]
    for line in eachline(stream)
        push!(nodes, [parse(T, e) for e in split(line)[1:3]])
    end
    nodes
end

#=
    Reads the triangles of a MSMS-generated surface mesh from the given .face file.

    @param stream
        Handle to .face file
    @param nodes
        List of reference nodes
    @param _
        Data type T for return value
    @return Vector{Triangle{T}}
=#
function readmsms_elements{T <: AbstractFloat}(stream::IOStream, nodes::Vector{Vector{T}}, ::Type{T}=Float64)
    elements = Triangle{T}[]
    # skip header lines
    [readline(stream) for i in 1:3]
    for line in eachline(stream)
        push!(elements, Triangle([nodes[parse(Int, e)] for e in split(line)[1:3]]...))
    end
    elements
end

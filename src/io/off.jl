#=
    Reads all data from the given OFF file.

    @param stream
        Handle to OFF file
    @param _
        Data type T for return value
    @return (Vector{Vector{T}}, Vector{Triangle{T}})
=#
function readoff{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    eof(stream) && return (Vector{T}[], Triangle{T}[])
    @assert readline(stream) == "OFF" "Invalid OFF file"

    # read number of nodes and elements
    numnodes, numelem = [parse(Int, s) for s in split(readline(stream))]

    nodes = readoff_nodes(stream, numnodes, T)
    (nodes, readoff_elements(stream, numelem, nodes, T))
end
readoff{T}(fname::String, ::Type{T}=Float64) = open(fh -> readoff(fh, T), fname)

#=
    Reads all node data from the given OFF file.

    @param stream
        Handle to OFF file
    @param n
        Total number of nodes
    @param _
        Data type T for return value
    @return Vector{Vector{T}}
=#
readoff_nodes{T <: AbstractFloat}(stream::IOStream, n::Int, ::Type{T}=Float64) = Vector{T}[[parse(T, e) for e in split(readline(stream))] for _ in 1:n]

#=
    Reads all element data from the given OFF file.

    @param stream
        Handle to OFF file
    @param n
        Total number of elements
    @param nodes
        List of reference nodes
    @param _
        Data type T for return value
    @return Vector{Triangle{T}}
=#
readoff_elements{T <: AbstractFloat}(stream::IOStream, n::Int, nodes::Vector{Vector{T}}, ::Type{T}=Float64) = Triangle{T}[Triangle([nodes[parse(Int, e) + 1] for e in split(readline(stream))[2:end]]...) for _ in 1:n]

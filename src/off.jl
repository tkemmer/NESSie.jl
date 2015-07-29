#=
    Reads all data from the given OFF file.

    @param stream
                Handle to OFF file
    @param dtype
                Data type T for return value
    @return (Vector{Vector{T}}, Vector{Element{T}})
=#
function readoff(stream::IOStream; dtype::Union(Type{Float64},Type{Float32})=Float64)
    @assert readline(stream) == "OFF\n" "Invalid OFF file"

    # read number of nodes and elements
    numnodes, numelem = [parse(Int, s) for s in split(readline(stream))]

    nodes = readoff_nodes(stream, numnodes, dtype)
    (nodes, readoff_elements(stream, numelem, nodes, dtype))
end

#=
    Reads all node data from the given OFF file.

    @param stream
                Handle to OFF file
    @param n
                Total number of nodes
    @param dtype
                Data type T for return value
    @return Vector{Vector{T}}
=#
function readoff_nodes(stream::IOStream, n::Int, dtype::Union(Type{Float64},Type{Float32})=Float64)
    nodes = Vector{dtype}[]

    @inbounds for _ in 1:n
        push!(nodes, [parse(dtype, e) for e in split(readline(stream))])
    end

    nodes
end


#=
    Reads all element data from the given OFF file.

    @param stream
                Handle to OFF file
    @param n
                Total number of elements
    @param nodes
                List of reference nodes
    @param dtype
                Data type T for return value
    @return Vector{Element{T}}
=#
function readoff_elements{T <: Union(Float64,Float32)}(stream::IOStream, n::Int, nodes::Vector{Vector{T}}, dtype::Union(Type{Float64},Type{Float32})=Float64)
    elements = Element{dtype}[]

    @inbounds for _ in 1:n
        push!(elements, Element([nodes[parse(Int, e) + 1] for e in split(readline(stream))[2:end]]...))
    end

    elements
end
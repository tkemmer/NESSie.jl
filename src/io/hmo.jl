#=
    Reads all data from the given HMO file.

    @param stream
                Handle to HMO file
    @param _
                Data type T for return value
    @return (Vector{Vector{T}}, Vector{Triangle{T}}, Vector{Charge{T}})
=#
function readhmo{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    nodes = readhmo_nodes(stream, T)
    (nodes, readhmo_elements(stream, nodes, T), readhmo_charges(stream, T))
end

#=
    Reads all node data from the given HMO file.

    @param stream
                Handle to HMO file
    @param _
                Data type T for return value
    @return Vector{Vector{T}}
=#
function readhmo_nodes{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    nodes = Vector{T}[]
    seek(stream, "BEG_NODL_DATA")
    readline(stream) # skip first line
    for line in eachline(stream)
        line == "END_NODL_DATA\n" && break
        push!(nodes, [parse(T, a) for a in split(line)[2:end]])
    end
    nodes
end

#=
    Reads all element data from the given HMO file.

    @param stream
                Handle to HMO file
    @param nodes
                List of referenced nodes
    @param _
                Data type T for return value
    @return Vector{Vector{T}}
=#
function readhmo_elements{T <: AbstractFloat}(stream::IOStream, nodes::Vector{Vector{T}}, ::Type{T}=T)
    elements = Triangle{T}[]
    seek(stream, "BEG_ELEM_DATA")
    readline(stream) # skip first line
    for line in eachline(stream)
        line == "END_ELEM_DATA\n" && break
        push!(elements, Triangle([nodes[parse(Int, a)] for a in split(line)[4:end]]...))
    end
    elements
end

#=
    Reads all charge data from the given HMO file.

    @param stream
                Handle to HMO file
    @param _
                Data type T for return value
    @return Vector{Charge{T}}
=#
function readhmo_charges{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    charges = Charge{T}[]
    seek(stream, "BEG_CHARGE_DATA")
    readline(stream) # skip first line
    for line in eachline(stream)
        line == "END_CHARGE_DATA\n" && break
        push!(charges, Charge([parse(T, a) for a in split(line)[2:end]]...))
    end
    charges
end

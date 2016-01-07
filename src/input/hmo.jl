#=
    Reads all data from the given HMO file.

    @param stream
                Handle to HMO file
    @param _
                Data type T for return value
    @return (Vector{Vector{T}}, Vector{Triangle{T}}, Vector{Charge{T}})
=#
function readhmo{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    nodes = Vector{T}[]
    elements = Triangle{T}[]
    charges = Charge{T}[]
    while !eof(stream)
        line = readline(stream)
        if line == "BEG_NODL_DATA\n"
            readline(stream)
            nodes = readhmo_nodes(stream, true, T)
        elseif line == "BEG_ELEM_DATA\n"
            readline(stream)
            elements = readhmo_elements(stream, nodes, true, T)
        elseif line == "BEG_CHARGE_DATA\n"
            readline(stream)
            charges = readhmo_charges(stream, true, T)
        end
    end
    (nodes, elements, charges)
end

#=
    Reads all node data from the given HMO file.

    @param stream
                Handle to HMO file
    @param atstart
                Specifies whether the handle points already to the first node
    @param _
                Data type T for return value
    @return Vector{Vector{T}}
=#
function readhmo_nodes{T <: AbstractFloat}(stream::IOStream, atstart::Bool=false, ::Type{T}=Float64)
    nodes = Vector{T}[]
    while !eof(stream)
        line = readline(stream)

        # If not already at the starting position, skip all lines
        # until we are.
        if !atstart
            if line == "BEG_NODL_DATA\n"
                atstart = true
                readline(stream)
            end
            continue
        end

        # Stop as soon as all nodes have been processed.
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
    @param atstart
                Specifies whether the handle points already to the first element
    @param _
                Data type T for return value
    @return Vector{Vector{T}}
=#
function readhmo_elements{T <: AbstractFloat}(stream::IOStream, nodes::Vector{Vector{T}}, atstart::Bool=false, ::Type{T}=T)
    elements = Triangle{T}[]
    while !eof(stream)
        line = readline(stream)

        # If not already at the starting position, skip all lines
        # until we are.
        if !atstart
            if line == "BEG_ELEM_DATA\n"
                atstart = true
                readline(stream)
            end
            continue
        end

        # Stop as soon as all elements have been processed.
        line == "END_ELEM_DATA\n" && break

        push!(elements, Triangle([nodes[parse(Int, a)] for a in split(line)[4:end]]...))
    end
    elements
end

#=
    Reads all charge data from the given HMO file.

    @param stream
                Handle to HMO file
    @param atstart
                Specifies whether the handle points already to the first charge
    @param _
                Data type T for return value
    @return Vector{Charge{T}}
=#
function readhmo_charges{T <: AbstractFloat}(stream::IOStream, atstart::Bool=false, ::Type{T}=Float64)
    charges = Charge{T}[]
    while !eof(stream)
        line = readline(stream)

        # If not already at the starting position, skip all lines
        # until we are.
        if !atstart
            if line == "BEG_CHARGE_DATA\n"
                atstart = true
                readline(stream)
            end
            continue
        end

        # Stop as soon as all charges have been processed.
        line == "END_CHARGE_DATA\n" && break

        push!(charges, Charge([parse(T, a) for a in split(line)[2:end]]...))
    end
    charges
end

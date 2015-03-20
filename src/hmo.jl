#=
    Reads all data from the given HMO file.

    @param stream
                Handle to HMO file
    @param dtype
                Data type T for return value
    @return (Vector{Vector{T}}, Vector{Element{T}}, Vector{Charge{T}})
=#
function readhmo(stream::IOStream; dtype::Union(Type{Float64},Type{Float32})=Float64)
    nodes = Vector{dtype}[]
    elements = Element{dtype}[]
    charges = Charge{dtype}[]
    while !eof(stream)
        line = readline(stream)
        if line == "BEG_NODL_DATA\n"
            readline(stream)
            nodes = readhmo_nodes(stream, true, dtype=dtype)
        elseif line == "BEG_ELEM_DATA\n"
            readline(stream)
            elements = readhmo_elements(stream, nodes, true, dtype=dtype)
        elseif line == "BEG_CHARGE_DATA\n"
            readline(stream)
            charges = readhmo_charges(stream, true, dtype=dtype)
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
    @param dtype
                Data type T for return value
    @return Vector{Vector{T}}
=#
function readhmo_nodes(stream::IOStream, atstart::Bool=false; dtype::Union(Type{Float64},Type{Float32})=Float64)
    nodes = Vector{dtype}[]
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
        if line == "END_NODL_DATA\n"
            break
        end

        push!(nodes, [parse(dtype, a) for a in split(line)[2:end]])
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
    @param dtype
                Data type T for return value
    @return Vector{Vector{T}}
=#
function readhmo_elements{T <: Union(Float64,Float32)}(stream::IOStream, nodes::Vector{Vector{T}}, atstart::Bool=false; dtype::Union(Type{Float64},Type{Float32})=T)
    elements = Element{dtype}[]
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
        if line == "END_ELEM_DATA\n"
            break
        end

        push!(elements, Element([nodes[parse(Int, a)] for a in split(line)[4:end]]...))
    end
    elements
end

#=
    Reads all charge data from the given HMO file.

    @param stream
                Handle to HMO file
    @param atstart
                Specifies whether the handle points already to the first charge
    @param dtype
                Data type T for return value
    @return Vector{Charge{T}}
=#
function readhmo_charges(stream::IOStream, atstart::Bool=false; dtype::Union(Type{Float64},Type{Float32})=Float64)
    charges = Charge{dtype}[]
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
        if line == "END_CHARGE_DATA\n"
            break
        end

        push!(charges, Charge(dtype, [parse(dtype, a) for a in split(line)[2:end]]...))
    end
    charges
end

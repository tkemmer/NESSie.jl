#=
    Reads protein partial charges from the given PQR file.
=#
function readpqr{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    charges = Charge{T}[]
    for line in eachline(stream)
        startswith(line, "ATOM") || continue # skip remarks and water molecules
        push!(charges, Charge([parse(T, e) for e in split(line)[end-4:end-1]]...))
    end
    charges
end

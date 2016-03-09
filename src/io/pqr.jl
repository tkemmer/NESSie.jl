#=
    Reads protein partial charges from the given PQR file.

    @param stream
        Handle to PQR file
    @param _
        Data type T for return value
=#
function readpqr{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    charges = Charge{T}[]
    for line in eachline(stream)
        startswith(line, "ATOM") || continue # skip remarks and water molecules
        push!(charges, Charge([parse(T, e) for e in split(line)[end-4:end-1]]...))
    end
    charges
end
readpqr{T}(fname::ASCIIString, ::Type{T}=Float64) = open(fh -> readpqr(fh, T), fname)

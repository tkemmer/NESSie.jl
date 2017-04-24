# =========================================================================================
"""
    readpqr{T <: AbstractFloat}(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads a charge model from the given PQR file.

# Return type
`Vector{Charge{T}}`

# Alias

    readpqr{T}(fname::String, ::Type{T}=Float64)

Reads the charge model using a file name rather than a `IOStream` object.
"""
function readpqr{T <: AbstractFloat}(stream::IOStream, ::Type{T}=Float64)
    charges = Charge{T}[]
    for line in eachline(stream)
        startswith(line, "ATOM") || continue # skip remarks and water molecules
        push!(charges, Charge([parse(T, e) for e in split(line)[end-4:end-1]]...))
    end
    charges
end
readpqr{T}(fname::String, ::Type{T}=Float64) = open(fh -> readpqr(fh, T), fname)

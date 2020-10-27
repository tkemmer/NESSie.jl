# =========================================================================================
"""
    readpqr{T <: AbstractFloat}(
        stream::IOStream,
              ::Type{T}=Float64
    )

Reads a charge model from the given PQR file.

# Specification
<https://pdb2pqr.readthedocs.io/en/latest/formats/pqr.html>

# Return type
`Vector{Charge{T}}`

# Alias

    readpqr{T}(fname::String, ::Type{T}=Float64)

Reads the charge model using a file name rather than a `IOStream` object.
"""
function readpqr(
        stream::IOStream,
              ::Type{T}=Float64
    ) where T <: AbstractFloat
    charges = Charge{T}[]
    for line in eachline(stream)
        startswith(line, "ATOM") || continue # skip remarks and water molecules
        vals = map(e -> parse(T, e), split(line)[end-4:end-1])
        vals[4] == 0 && continue             # skip zero charges
        push!(charges, Charge(vals...))
    end
    charges
end

function readpqr(
        fname::String,
             ::Type{T}=Float64
    ) where T
    open(fh -> readpqr(fh, T), fname)
end

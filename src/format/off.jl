# =========================================================================================
"""
    readoff(stream::IOStream, ::Type{T}=Float64)

Reads a surface model from the given OFF file.

!!! note
    This file type does not support charge models! Hence, the charge list of the returning
    `Model` object is empty and has to be set separately.

# Specification
<http://www.geomview.org/docs/html/OFF.html>

# Return type
[`Model{T}`](@ref)

# Alias

    readoff(fname::AbstractString, ::Type{T}=Float64)

Reads the model using a file name rather than a `IOStream` object.
"""
@inline function readoff(stream::IOStream, ::Type{T}=Float64) where T
    _readoff(Stream{format"OFF"}(stream), T)
end

@inline function readoff(fname::AbstractString, ::Type{T}=Float64) where T
    _readoff(File{format"OFF"}(fname), T)
end

@inline function _readoff(f, ::Type{T}) where T <: AbstractFloat
    Model(load(f; pointtype = _pointtype(T)))
end

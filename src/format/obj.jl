# =========================================================================================
"""
    readobj(fname::AbstractString, ::Type{T}=Float64)

Reads a surface model from the given Wavefront OBJ file.

!!! note
    This file type does not support charge models! Hence, the charge list of the returning
    `Model` object is empty and has to be set separately.

# Specification
- <https://paulbourke.net/dataformats/obj/>

# Return type
[`Model{T, Triangle{T}}`](@ref)
"""
@inline function readobj(fname::AbstractString, ::Type{T}=Float64) where T
    _readobj(File{format"OBJ"}(fname), T)
end

@inline function _readobj(f, ::Type{T}) where T <: AbstractFloat
    Model(load(f; pointtype = _pointtype(T)))
end


# =========================================================================================
"""
    function writeobj(stream::IOStream, model ::Model{T, Triangle{T}})

Creates a Wavefront OBJ file from a given surface model.

# Specification
- <https://paulbourke.net/dataformats/obj/>

# Return type
`Nothing`

# Alias

    writeobj(fname::AbstractString, model::Model{T, Triangle{T}})

Creates the OBJ file by name rather than `IOStream` object.
"""
@inline function writeobj(stream::IOStream, model::Model{T, Triangle{T}}) where T
    _writeobj(Stream{format"OBJ"}(stream), model)
end

@inline function writeobj(fname::AbstractString, model::Model{T, Triangle{T}}) where T
    _writeobj(File{format"OBJ"}(fname), model)
end

@inline function _writeobj(f, model::Model{T, Triangle{T}}) where T
    save(f, GeometryBasics.normal_mesh(GeometryBasics.mesh(model)))
end

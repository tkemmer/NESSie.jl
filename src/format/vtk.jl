# =========================================================================================
"""
    writevtk(stream::IOStream, model::Model{T, Triangle{T}})

Creates a VTK PolyData file from a given surface Model.

# Specification
<https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf>

# Return type
`Nothing`

# Alias

    writevtk(fname::AbstractString, model::Model{T, Triangle{T}})

Creates the VTK file by name rather than `IOStream` object. The file extension (`.vtp`) is
appended automatically if not present.
"""
function writevtk(
        stream::IOStream,
        model ::Model{T, Triangle{T}}
    ) where T
    mktempdir() do dir
        fname = joinpath(dir, "mesh")
        outfiles = vtk_save(_vtk_grid(fname, model))
        open(only(outfiles), "r") do io
            write(stream, io)
        end
    end
    nothing
end

@inline function writevtk(
        fname::AbstractString,
        model::Model{T, Triangle{T}}
    ) where T
    vtk_save(_vtk_grid(fname, model))
    nothing
end

@inline function _vtk_grid(
        fname::AbstractString,
        model ::Model{T, Triangle{T}}
    ) where T
    revidx = _reverseindex(model.nodes)

    points = Matrix{T}(undef, 3, length(model.nodes))
    for (i, node) in enumerate(model.nodes)
        points[:, i] = node
    end

    cells = Vector{MeshCell{PolyData.Polys, Vector{Int32}}}(undef, length(model.elements))
    for (i, elem) in enumerate(model.elements)
        cells[i] = MeshCell(
            PolyData.Polys(),
            Int32[revidx[elem.v1], revidx[elem.v2], revidx[elem.v3]]
        )
    end

    vtk_grid(fname, points, cells; append = false, ascii = true)
end

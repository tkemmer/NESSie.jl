# =========================================================================================
"""
    writevtk(fname::AbstractString, model::Model{T, Triangle{T}})

Creates a VTK PolyData file from a given surface Model. The file extension (`.vtp`) is
appended automatically if not present.

# Specification
<https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf>

# Supported keyword arguments
 - `ascii::Bool = true`: use `ascii` instead of `binary` as data array format

# Return type
`Nothing`
"""
@inline function writevtk(
    fname::AbstractString,
    model::Model{T, Triangle{T}};
    ascii::Bool = true
) where T
    vtk_save(_vtk_grid(fname, model; ascii))
    nothing
end

@inline function _vtk_grid(
    fname::AbstractString,
    model ::Model{T, Triangle{T}};
    ascii::Bool = true
) where T
    revidx = _reverseindex(model.nodes)

    points = Matrix{T}(undef, 3, length(model.nodes))
    for (i, node) in enumerate(model.nodes)
        points[:, i] = node
    end

    cells = [
        MeshCell(
            PolyData.Polys(),
            Int32[revidx[elem.v1], revidx[elem.v2], revidx[elem.v3]]
        ) for elem in model.elements
    ]

    vtk_grid(fname, points, cells; append = false, ascii)
end

# =========================================================================================
"""
    writevtk(
        stream::IOStream,
        model ::Model{T}
    )

Creates a VTK-compatible output file from a given surface or volume model. The exact file
type is determined by the given model:

| Model type    | Resulting file type  |
|---------------|----------------------|
| Surface model | VTK PolyData         |
| Volume model  | VTK UnstructuredGrid |

# Specification
<https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf>

# Return type
`Void`

# Alias

    writevtk(
        fname::String,
        model::Model{T}
    )

Creates the VTK file by name rather than `IOStream` object.
"""
function writevtk(
        stream::IOStream,
        model ::Model{T, Triangle{T}}
    ) where T
    xdoc = XMLDocument()
    revidx = _reverseindex(model.nodes)

    # VTKFile
    xroot = create_root(xdoc, "VTKFile")
    set_attribute(xroot, "type", "PolyData")

    # PolyData/Piece
    xpiece = new_child(new_child(xroot, "PolyData"), "Piece")
    set_attribute(xpiece, "NumberOfPoints", length(model.nodes))
    set_attribute(xpiece, "NumberOfVerts",  0)
    set_attribute(xpiece, "NumberOfLines",  0)
    set_attribute(xpiece, "NumberOfStrips", 0)
    set_attribute(xpiece, "NumberOfPolys", length(model.elements))

    # Points
    xpoints = new_child(new_child(xpiece, "Points"), "DataArray")
    set_attribute(xpoints, "type", "$T")
    set_attribute(xpoints, "NumberOfComponents", "3")
    set_attribute(xpoints, "format", "ascii")
    add_text(xpoints, join(Iterators.flatten(model.nodes), " "))

    # Polys
    xpolys = new_child(xpiece, "Polys")

    # Polys/connectivity
    xconn = new_child(xpolys, "DataArray")
    set_attribute(xconn, "type", "Int32")
    set_attribute(xconn, "Name", "connectivity")
    set_attribute(xconn, "format", "ascii")
    add_text(xconn, join([revidx[n]-1 for n in
        Iterators.flatten(Vector{T}[o.v1, o.v2, o.v3] for o in model.elements)], " "))

    # Polys/offsets
    xoffs = new_child(xpolys, "DataArray")
    set_attribute(xoffs, "type", "Int32")
    set_attribute(xoffs, "Name", "offsets")
    set_attribute(xoffs, "format", "ascii")
    add_text(xoffs, join([3i for i in 1:length(model.elements)], " "))

    # create .vtp file
    println(stream, string(xdoc))
    nothing
end

function writevtk(
        stream::IOStream,
        model ::Model{T, Tetrahedron{T}}
    ) where T
    xdoc = XMLDocument()
    revidx = _reverseindex(model.nodes)

    # VTKFile
    xroot = create_root(xdoc, "VTKFile")
    set_attribute(xroot, "type", "UnstructuredGrid")

    # UnstructuredGrid/Piece
    xpiece = new_child(new_child(xroot, "UnstructuredGrid"), "Piece")
    set_attribute(xpiece, "NumberOfPoints", length(model.nodes))
    set_attribute(xpiece, "NumberOfCells",  length(model.elements))

    # Points
    xpoints = new_child(new_child(xpiece, "Points"), "DataArray")
    set_attribute(xpoints, "type", "$T")
    set_attribute(xpoints, "NumberOfComponents", "3")
    set_attribute(xpoints, "format", "ascii")
    add_text(xpoints, join(Iterators.flatten(model.nodes), " "))

    # Cells
    xcells = new_child(xpiece, "Cells")

    # Cells/connectivity
    xconn = new_child(xcells, "DataArray")
    set_attribute(xconn, "type", "Int32")
    set_attribute(xconn, "Name", "connectivity")
    set_attribute(xconn, "format", "ascii")
    add_text(xconn, join([revidx[n]-1 for n in
        Iterators.flatten(Vector{T}[o.v1, o.v2, o.v3, o.v4] for o in model.elements)], " "))

    # Cells/offsets
    xoffs = new_child(xcells, "DataArray")
    set_attribute(xoffs, "type", "Int32")
    set_attribute(xoffs, "Name", "offsets")
    set_attribute(xoffs, "format", "ascii")
    add_text(xoffs, join([4i for i in 1:length(model.elements)], " "))

    # Cells/types
    xtypes = new_child(xcells, "DataArray")
    set_attribute(xtypes, "type", "UInt8")
    set_attribute(xtypes, "Name", "types")
    set_attribute(xtypes, "format", "ascii")
    add_text(xtypes, "10 "^length(model.elements))

    # create .vtu file
    println(stream, string(xdoc))
    nothing
end

@inline function writevtk(
        fname::String,
        model::M
    ) where {T, M <: Model{T}}
    open(fh -> writevtk(fh, model), fname, "w")
end

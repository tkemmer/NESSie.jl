#=
    Creates a VTK UnstructuredGrid file from a given volume model.

    File format specification:
    http://www.vtk.org/VTK/img/file-formats.pdf

    @param fname/stream
        Path or handle to (writable) VTK file
    @param model
        A volume model
    @return nothing
=#
function writevtk{T}(stream::IOStream, model::VolumeModel{T})
    xdoc = XMLDocument()
    revidx = reverseindex(model.nodes)

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
    add_text(xpoints, join(unpack(model.nodes), " "))

    # Cells
    xcells = new_child(xpiece, "Cells")

    # Cells/connectivity
    xconn = new_child(xcells, "DataArray")
    set_attribute(xconn, "type", "Int32")
    set_attribute(xconn, "Name", "connectivity")
    set_attribute(xconn, "format", "ascii")
    add_text(xconn, join([revidx[object_id(n)]-1 for n in unpack([Vector{T}[o.v1, o.v2, o.v3, o.v4] for o in model.elements])], " "))

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
writevtk{T}(fname::String, model::VolumeModel{T}) = open(fh -> writevtk(fh, model), fname, "w")

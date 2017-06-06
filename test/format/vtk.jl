using NESSie.Format
using LightXML: parse_string, root, name, child_elements, attribute, content

context("writevtk (surface model)") do
    for T in testtypes
        # empty system
        model = Model{T, Triangle{T}}()
        (fname, fh) = mktemp()
        try
            xroot = root(parse_string(readback(fh -> writevtk(fh, model))))
            @fact name(xroot) --> "VTKFile"
            @fact attribute(xroot, "type") --> "PolyData"
            xgrid = collect(child_elements(xroot))
            @fact length(xgrid) --> 1
            @fact name(xgrid[1]) --> "PolyData"
            xgrid = collect(child_elements(xgrid[1]))
            @fact length(xgrid) --> 1
            @fact name(xgrid[1]) --> "Piece"
            @fact attribute(xgrid[1], "NumberOfPoints") --> "0"
            @fact attribute(xgrid[1], "NumberOfVerts") --> "0"
            @fact attribute(xgrid[1], "NumberOfLines") --> "0"
            @fact attribute(xgrid[1], "NumberOfStrips") --> "0"
            @fact attribute(xgrid[1], "NumberOfPolys") --> "0"
            xpiece = collect(child_elements(xgrid[1]))
            @fact length(xpiece) --> 2
            @fact name(xpiece[1]) --> "Points"
            xdata = collect(child_elements(xpiece[1]))
            @fact length(xdata) --> 1
            @fact name(xdata[1]) --> "DataArray"
            @fact attribute(xdata[1], "type") --> "$T"
            @fact attribute(xdata[1], "NumberOfComponents") --> "3"
            @fact attribute(xdata[1], "format") --> "ascii"
            @fact strip(content(xdata[1])) --> ""
            @fact name(xpiece[2]) --> "Polys"
            xdata = collect(child_elements(xpiece[2]))
            @fact length(xdata) --> 2
            @fact name(xdata[1]) --> "DataArray"
            @fact attribute(xdata[1], "type") --> "Int32"
            @fact attribute(xdata[1], "Name") --> "connectivity"
            @fact attribute(xdata[1], "format") --> "ascii"
            @fact strip(content(xdata[1])) --> ""
            @fact name(xdata[2]) --> "DataArray"
            @fact attribute(xdata[2], "type") --> "Int32"
            @fact attribute(xdata[2], "Name") --> "offsets"
            @fact attribute(xdata[2], "format") --> "ascii"
            @fact strip(content(xdata[2])) --> ""
        finally
            close(fh)
            rm(fname)
        end

        # small system
        nodes = Vector{T}[[0, 0, 0], [0, 0, 3], [0, 3, 0], [1, -3, 3]]
        model = Model(
            nodes,
            Triangle{T}[
                Triangle(nodes[1], nodes[2], nodes[3]),
                Triangle(nodes[1], nodes[4], nodes[2])
            ]
        )
        (fname, fh) = mktemp()
        try
            xroot = root(parse_string(readback(fh -> writevtk(fh, model))))
            @fact name(xroot) --> "VTKFile"
            @fact attribute(xroot, "type") --> "PolyData"
            xgrid = collect(child_elements(xroot))
            @fact length(xgrid) --> 1
            @fact name(xgrid[1]) --> "PolyData"
            xgrid = collect(child_elements(xgrid[1]))
            @fact length(xgrid) --> 1
            @fact name(xgrid[1]) --> "Piece"
            @fact attribute(xgrid[1], "NumberOfPoints") --> "4"
            @fact attribute(xgrid[1], "NumberOfVerts") --> "0"
            @fact attribute(xgrid[1], "NumberOfLines") --> "0"
            @fact attribute(xgrid[1], "NumberOfStrips") --> "0"
            @fact attribute(xgrid[1], "NumberOfPolys") --> "2"
            xpiece = collect(child_elements(xgrid[1]))
            @fact length(xpiece) --> 2
            @fact name(xpiece[1]) --> "Points"
            xdata = collect(child_elements(xpiece[1]))
            @fact length(xdata) --> 1
            @fact name(xdata[1]) --> "DataArray"
            @fact attribute(xdata[1], "type") --> "$T"
            @fact attribute(xdata[1], "NumberOfComponents") --> "3"
            @fact attribute(xdata[1], "format") --> "ascii"
            @fact [parse(T, e) for e in split(strip(content(xdata[1])))] --> T[0, 0, 0, 0, 0, 3, 0, 3, 0, 1, -3, 3]
            @fact name(xpiece[2]) --> "Polys"
            xdata = collect(child_elements(xpiece[2]))
            @fact length(xdata) --> 2
            @fact name(xdata[1]) --> "DataArray"
            @fact attribute(xdata[1], "type") --> "Int32"
            @fact attribute(xdata[1], "Name") --> "connectivity"
            @fact attribute(xdata[1], "format") --> "ascii"
            @fact [parse(Int, e) for e in split(strip(content(xdata[1])))] --> [0, 1, 2, 0, 3, 1]
            @fact name(xdata[2]) --> "DataArray"
            @fact attribute(xdata[2], "type") --> "Int32"
            @fact attribute(xdata[2], "Name") --> "offsets"
            @fact attribute(xdata[2], "format") --> "ascii"
            @fact [parse(Int, e) for e in split(strip(content(xdata[2])))] --> [3, 6]
        finally
            close(fh)
            rm(fname)
        end
    end
end

context("writevtk (volume model)") do
    for T in testtypes
        # empty system
        model = Model{T, Tetrahedron{T}}()
        (fname, fh) = mktemp()
        try
            xroot = root(parse_string(readback(fh -> writevtk(fh, model))))
            @fact name(xroot) --> "VTKFile"
            @fact attribute(xroot, "type") --> "UnstructuredGrid"
            xgrid = collect(child_elements(xroot))
            @fact length(xgrid) --> 1
            @fact name(xgrid[1]) --> "UnstructuredGrid"
            xgrid = collect(child_elements(xgrid[1]))
            @fact length(xgrid) --> 1
            @fact name(xgrid[1]) --> "Piece"
            @fact attribute(xgrid[1], "NumberOfPoints") --> "0"
            @fact attribute(xgrid[1], "NumberOfCells") --> "0"
            xpiece = collect(child_elements(xgrid[1]))
            @fact length(xpiece) --> 2
            @fact name(xpiece[1]) --> "Points"
            xdata = collect(child_elements(xpiece[1]))
            @fact length(xdata) --> 1
            @fact name(xdata[1]) --> "DataArray"
            @fact attribute(xdata[1], "type") --> "$T"
            @fact attribute(xdata[1], "NumberOfComponents") --> "3"
            @fact attribute(xdata[1], "format") --> "ascii"
            @fact strip(content(xdata[1])) --> ""
            @fact name(xpiece[2]) --> "Cells"
            xdata = collect(child_elements(xpiece[2]))
            @fact length(xdata) --> 3
            @fact name(xdata[1]) --> "DataArray"
            @fact attribute(xdata[1], "type") --> "Int32"
            @fact attribute(xdata[1], "Name") --> "connectivity"
            @fact attribute(xdata[1], "format") --> "ascii"
            @fact strip(content(xdata[1])) --> ""
            @fact name(xdata[2]) --> "DataArray"
            @fact attribute(xdata[2], "type") --> "Int32"
            @fact attribute(xdata[2], "Name") --> "offsets"
            @fact attribute(xdata[2], "format") --> "ascii"
            @fact strip(content(xdata[2])) --> ""
            @fact name(xdata[3]) --> "DataArray"
            @fact attribute(xdata[3], "type") --> "UInt8"
            @fact attribute(xdata[3], "Name") --> "types"
            @fact attribute(xdata[3], "format") --> "ascii"
            @fact strip(content(xdata[3])) --> ""
        finally
            close(fh)
            rm(fname)
        end

        # small system
        nodes = Vector{T}[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [-2, 0, 0]]
        model = Model(
            nodes,
            Tetrahedron{T}[
                Tetrahedron(nodes[1], nodes[2], nodes[3], nodes[4]),
                Tetrahedron(nodes[1], nodes[2], nodes[3], nodes[5])
            ]
        )
        (fname, fh) = mktemp()
        try
            xroot = root(parse_string(readback(fh -> writevtk(fh, model))))
            @fact name(xroot) --> "VTKFile"
            @fact attribute(xroot, "type") --> "UnstructuredGrid"
            xgrid = collect(child_elements(xroot))
            @fact length(xgrid) --> 1
            @fact name(xgrid[1]) --> "UnstructuredGrid"
            xgrid = collect(child_elements(xgrid[1]))
            @fact length(xgrid) --> 1
            @fact name(xgrid[1]) --> "Piece"
            @fact attribute(xgrid[1], "NumberOfPoints") --> "5"
            @fact attribute(xgrid[1], "NumberOfCells") --> "2"
            xpiece = collect(child_elements(xgrid[1]))
            @fact length(xpiece) --> 2
            @fact name(xpiece[1]) --> "Points"
            xdata = collect(child_elements(xpiece[1]))
            @fact length(xdata) --> 1
            @fact name(xdata[1]) --> "DataArray"
            @fact attribute(xdata[1], "type") --> "$T"
            @fact attribute(xdata[1], "NumberOfComponents") --> "3"
            @fact attribute(xdata[1], "format") --> "ascii"
            @fact [parse(T, e) for e in split(strip(content(xdata[1])))] --> T[0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, -2, 0, 0]
            @fact name(xpiece[2]) --> "Cells"
            xdata = collect(child_elements(xpiece[2]))
            @fact length(xdata) --> 3
            @fact name(xdata[1]) --> "DataArray"
            @fact attribute(xdata[1], "type") --> "Int32"
            @fact attribute(xdata[1], "Name") --> "connectivity"
            @fact attribute(xdata[1], "format") --> "ascii"
            @fact [parse(Int, e) for e in split(strip(content(xdata[1])))] --> [0, 1, 2, 3, 0, 1, 2, 4]
            @fact name(xdata[2]) --> "DataArray"
            @fact attribute(xdata[2], "type") --> "Int32"
            @fact attribute(xdata[2], "Name") --> "offsets"
            @fact attribute(xdata[2], "format") --> "ascii"
            @fact [parse(Int, e) for e in split(strip(content(xdata[2])))] --> [4, 8]
            @fact name(xdata[3]) --> "DataArray"
            @fact attribute(xdata[3], "type") --> "UInt8"
            @fact attribute(xdata[3], "Name") --> "types"
            @fact attribute(xdata[3], "format") --> "ascii"
            @fact [parse(Int, e) for e in split(strip(content(xdata[3])))] --> [10, 10]
        finally
            close(fh)
            rm(fname)
        end
    end
end
